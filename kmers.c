#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pthread.h>
#include<ctype.h>
#include<unistd.h>

#include"dataset.h"
#include"kmers.h"

#ifdef __SSE2__
#include<emmintrin.h>
#endif

typedef struct {
  kmer_frequencies* freq;
  char** kmers;
  dataset* ds;
  size_t* kmer_length;
  int thread_num;
  int* n_threads;
} thread_handle;

pthread_mutex_t lock;

static inline char* gen_dna_string_from_binary(size_t kmer_length,
					       unsigned long binary) {

  char* s = (char*)malloc(sizeof(char)*(kmer_length+1));

  unsigned long i;
  
  unsigned long test_a;
  unsigned long test_b;

  unsigned long res;

  for(i=0;i<kmer_length;i++) {

    res = 0;
    
    test_a = (1 << (i*2));
    test_b = (1 << (i*2)+1);

    if(test_a == (binary & test_a)) res += 1;
    if(test_b == (binary & test_b)) res += 2;

    switch (res) {
    case 0:
      s[i] = 'A';
      break;
    case 1:
      s[i] = 'C';
      break;
    case 2:
      s[i] = 'G';
      break;
    case 3:
      s[i] = 'T';
      break;
    }
  }
  s[kmer_length] = 0;
  return(s);
}

size_t number_of_kmers_from_length (size_t kmer_length) {

  int i;
  size_t n_kmer = 4;
  for(i = 0; i < kmer_length-1; i++) {
    n_kmer *= 4;
  }
  return (n_kmer);
}

size_t number_of_protein_kmers_from_length(size_t protein_kmer_length) {
  int i;
  size_t n_kmer = 20;
  for(i = 0; i < protein_kmer_length-1; i++) {
    n_kmer *= 20;
  }
  return (n_kmer);
}

char** gen_dna_kmers(size_t kmer_length) {

  unsigned long i;

  unsigned long a = 0;

  char** kmers;
  
  for(i=0;i<kmer_length*2;i++) {
    a |= (1 << i);
  }

  //  printf("generating %i kmers\n", a+1);

  kmers = (char**)malloc(sizeof(char*)*(a+1));
  
  for(i=a;i>0;i--) {
    kmers[i] = gen_dna_string_from_binary(kmer_length, i);
  }
  kmers[0] = gen_dna_string_from_binary(kmer_length, 0);
  return(kmers);
}

void recursive_protein_kmer_generator(char** kmers,
				      size_t length, size_t level,
				      char* kmer_template,
				      char* lookup_table,
				      int* global_counter) {
  int i;
  if (level < length) {
    level++;
    for(i=0;i<20;i++) {
      kmer_template[level-1] = lookup_table[i];
      recursive_protein_kmer_generator(kmers, length, level,
				       kmer_template,lookup_table,
				       global_counter);
    }
    level--;
  } else {
    kmer_template[level] = lookup_table[i];
    kmers[global_counter[0]] = (char*)malloc(sizeof(char)*length);
    memcpy(kmers[global_counter[0]], kmer_template, sizeof(char)*length);
    global_counter[0]++;
  }
}	
      
char** gen_protein_kmers(size_t protein_kmer_length) {

  char lookup_table[20] = {'A','C','D','E','F',
			   'G','H','I','K','L',
			   'M','N','P','Q','R',
			   'S','T','V','W','Y'};

  char* kmer_template = (char*)malloc(sizeof(char)*protein_kmer_length);

  int global_counter = 0;
  
  char** kmers =
    (char**)malloc(sizeof(char*)
		   *number_of_protein_kmers_from_length(protein_kmer_length));
  
  recursive_protein_kmer_generator(kmers,
				   protein_kmer_length,
				   0,
				   kmer_template,
				   lookup_table,
				   &global_counter);

  return(kmers);
}

short int frequence_of_kmer_in_sequence(char* kmer, char* sequence,
					size_t kmer_length) {

  int sequence_length = strlen(sequence);
  int comparisions;
  int i,j;

  int control;

  short int matches = 0;
  
  if ( kmer_length > sequence_length) return 0;

  comparisions = sequence_length-kmer_length;
  
#ifdef __SSE2__ 
  int v_comparisons = sequence_length - 16;
  unsigned char kmer_c_vector[16] __attribute__((aligned(16)));
  unsigned char kmer_c_mask[16] __attribute__((aligned(16)));

  __m128i kmer_vector;
  __m128i kmer_mask;
  __m128i current_sequence;
  __m128i comparison_result;
  int magic_return;
  int win = 65535; /* 2^16-1 */

  if ( (v_comparisons < 0) || (kmer_length > 16) ) {
    v_comparisons = 0;
  } else {
    for(i = 0; i<kmer_length;i++) {
      kmer_c_vector[i] = kmer[i];
      kmer_c_mask[i] = 255;
    }
    for(i= kmer_length ; i<16;i++) {
      kmer_c_vector[i] = 0;
      kmer_c_mask[i] = 0;
    }
    
    kmer_vector = _mm_load_si128(kmer_c_vector);
    kmer_mask = _mm_load_si128(kmer_c_mask);
    
    for(i=0; i<v_comparisons;i++) {
      current_sequence = _mm_loadu_si128((sequence+i));
      current_sequence = _mm_and_si128(current_sequence,kmer_mask);
      comparison_result = _mm_cmpeq_epi8(current_sequence,kmer_vector);
      magic_return = _mm_movemask_epi8(comparison_result);
      matches += (win == magic_return);
    }
  }
  /* handle not vectorizeable tail */
  for(i = v_comparisons; i<comparisions; i++) {
    control = 0;
    for(j = 0; j<kmer_length; j++) {
      if ( kmer[j] != sequence[j+i] ) {
	control = 1;
	goto out_of_inner_for;
      }
    }
  out_of_inner_for:
    if(control == 0) {
      matches++;
    }
  }
  
#else
  
  for(i = 0; i<comparisions; i++) {
    control = 0;
    for(j = 0; j<kmer_length; j++) {
      if ( kmer[j] != sequence[j+i] ) {
	control = 1;
	goto out_of_inner_for;
      }
    }
  out_of_inner_for:
    if(control == 0) {
      matches++;
    }
  }
  
#endif
   
  return (matches);
}

void free_kmers(char** kmers, size_t kmer_length, int protein) {

  size_t i;
  size_t n;
  if(protein) {
    n = number_of_protein_kmers_from_length(kmer_length);
  } else {
    n = number_of_kmers_from_length(kmer_length);
  }
  for(i=0;i<n;i++) {
    free(kmers[i]);
  }

  free(kmers);
}

void* frequencies_from_dataset_thread_handler(void* arg) {
  int i,j;

  pthread_mutex_lock(&lock);
  
  thread_handle* th = (thread_handle*) arg;
  size_t stride = (th->freq->n_seq)/(th->n_threads[0]);
  size_t start_seq = th->thread_num*stride;
  size_t end_seq = start_seq+stride;
  
  size_t kmer_length = th->kmer_length[0];

  short int** frequencies = th->freq->frequencies;
  char** sequences = th->ds->sequences;
  char** kmers = th->kmers;
  int n_kmers = th->freq->n_kmers;

  pthread_mutex_unlock(&lock);
  
  for(i=start_seq;i<end_seq;i++) {
    for(j=0;j<n_kmers;j++) {
      frequencies[i][j] =
	frequence_of_kmer_in_sequence(kmers[j],
				      sequences[i],
				      kmer_length);
    }
  }
}

kmer_frequencies frequencies_from_dataset(dataset ds, size_t kmer_length,
					  int n_threads, int protein) {

  int i,j;

  kmer_frequencies freq;

  char** kmers;
  char*** kmers_t;
  
  thread_handle* th = (thread_handle*)malloc(sizeof(thread_handle)*n_threads);
  pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);
  size_t stride;

  if(protein) {
    kmers = gen_protein_kmers(kmer_length);
    freq.n_kmers = number_of_protein_kmers_from_length(kmer_length);
  } else {
    kmers = gen_dna_kmers(kmer_length);
    freq.n_kmers = number_of_kmers_from_length(kmer_length);
  }
  freq.n_seq = ds.n_values;
  
  freq.frequencies = (short int**)malloc(sizeof(short int*)*freq.n_seq);

  for(i=0;i<freq.n_seq;i++) {
    freq.frequencies[i] = (short int*)malloc(sizeof(short int)*freq.n_kmers);
  }

  if( n_threads < freq.n_seq ) {

    pthread_mutex_init(&lock,NULL);
    
    kmers_t = (char***)malloc(sizeof(char**)*n_threads);
  
    for(i=0;i<n_threads;i++) {

      if(protein) {
	kmers_t[i] = gen_protein_kmers(kmer_length);
      } else {
	kmers_t[i] = gen_dna_kmers(kmer_length);
      }
      
      th[i].freq = &freq;
      th[i].kmers = kmers_t[i];
      th[i].ds = &ds;
      th[i].kmer_length = &kmer_length;
      th[i].thread_num = i;
      th[i].n_threads = &n_threads;

      pthread_create(threads+i, NULL,
		     frequencies_from_dataset_thread_handler,
		     th+i);
      
    }

    for(i=0;i<n_threads;i++) {
      pthread_join(threads[i], NULL);
      free_kmers(kmers_t[i],kmer_length,protein);
      pthread_mutex_destroy(&lock);
    }

    stride = freq.n_seq/n_threads;
    
    for(i=stride*n_threads;i<(stride*n_threads)+(freq.n_seq%n_threads);i++) {
      for(j=0;j<freq.n_kmers;j++) {
	freq.frequencies[i][j] =
	  frequence_of_kmer_in_sequence(kmers[j],ds.sequences[i],kmer_length);
      }
    }

    free(kmers_t);
    
  } else {

    for(i=0;i<freq.n_seq;i++) {
      for(j=0;j<freq.n_kmers;j++) {
	freq.frequencies[i][j] =
	  frequence_of_kmer_in_sequence(kmers[j],ds.sequences[i],kmer_length);
      }
    }
  }
  
  free_kmers(kmers,kmer_length,protein);
  free(th);
  free(threads);
  
  return(freq);
}

void write_kmer_base(FILE* f,kmer_frequencies freq) {

  int i,j;

  for(i=0;i<freq.n_seq;i++) {
    fprintf(f,"sequence_%i\t",i);
    for(j=0;j<freq.n_kmers-1;j++) {
      fprintf(f,"%hi\t",freq.frequencies[i][j]);
    }
    fprintf(f,"%hi\n",freq.frequencies[i][freq.n_kmers-1]);
  }
}

void free_kmer_frequencies(kmer_frequencies freq) {

  int i;
  for(i = 0; i<freq.n_seq;i++) {
    free(freq.frequencies[i]);
  }
  free(freq.frequencies);
}

dataset generate_n_permutations_for_sequence(dataset ds, size_t index) {

  dataset return_set;

  size_t length = ds.sequence_lengths[index];
  char* sequence = ds.sequences[index];

  size_t i,j;
  size_t N_count = 0, n_permutations;
  size_t* n_position_array = (size_t*)malloc(sizeof(size_t)*length);

  char** kmers;

  char** sequences;
  
  if(index > ds.n_values) {
    printf("Error index out of range in"
	   "generate_n_permultations_for_sequence()\n");
    _exit(1);
  }
  
  for(i=0;i<length;i++) {
    if ( 'N' == toupper(sequence[i]) ) {
      n_position_array[N_count] = i;
      N_count++;
    }
  }
  
  if(N_count > sizeof(unsigned long)*4) {
    printf("Error to many Ns in sequence \n"
	   "generate_n_permultations_for_sequence()\n"
	   "is currently limited to %i Ns in sequence\n",
	   (int)sizeof(unsigned long)*4);
    _exit(1);
  }
 
  kmers = gen_dna_kmers(N_count);
  
  n_permutations = number_of_kmers_from_length(N_count);
  sequences = (char**)malloc(sizeof(char*)*n_permutations);

  return_set.sequence_lengths = (size_t*)malloc(sizeof(size_t)*n_permutations);
  
  for(i=0;i<n_permutations;i++) {
    sequences[i] = (char*)malloc(sizeof(char)*length+1);
    memcpy(sequences[i],sequence,length);
    sequences[i][length]=0;
    for(j=0;j<N_count;j++) {
      sequences[i][n_position_array[j]] = kmers[i][j];
    }
    return_set.sequence_lengths[i] = length;
  }
  free_kmers(kmers,N_count,0);

  return_set.sequences = sequences;
  return_set.max_sequence_length = length;
  return_set.n_values = n_permutations;

  return(return_set);
 
}
		
  
  
