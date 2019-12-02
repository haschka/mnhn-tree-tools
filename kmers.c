#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pthread.h>

#include"dataset.h"
#include"kmers.h"

typedef struct {
  kmer_frequencies* freq;
  char** kmers;
  dataset* ds;
  size_t* kmer_length;
  int thread_num;
  int* n_threads;
} thread_handle;

pthread_mutex_t lock;

static inline char* gen_dna_string_from_binary(size_t kmer_length, int binary) {

  char* s = (char*)malloc(sizeof(char)*(kmer_length+1));

  int i;
  
  int test_a;
  int test_b;

  int res;

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

char** gen_dna_kmers(size_t kmer_length) {

  int i;

  int a = 0;

  char** kmers;
  
  for(i=0;i<kmer_length*2;i++) {
    a |= (1 << i);
  }

  //  printf("generating %i kmers\n", a+1);

  kmers = (char**)malloc(sizeof(char*)*(a+1));
  
  for(i=a;i>=0;i--) {
    kmers[i] = gen_dna_string_from_binary(kmer_length, i);
  }
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
  return (matches);
}

void free_kmers(char** kmers, size_t kmer_length) {

  int i;
  int n = number_of_kmers_from_length(kmer_length);

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
					  int n_threads) {

  int i,j;

  kmer_frequencies freq;

  char** kmers = gen_dna_kmers(kmer_length);

  char*** kmers_t;
  
  thread_handle* th = (thread_handle*)malloc(sizeof(thread_handle)*n_threads);
  pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);
  size_t stride;
  
  freq.n_kmers = number_of_kmers_from_length(kmer_length);
  freq.n_seq = ds.n_values;
  
  freq.frequencies = (short int**)malloc(sizeof(short int*)*freq.n_seq);

  for(i=0;i<freq.n_seq;i++) {
    freq.frequencies[i] = (short int*)malloc(sizeof(short int)*freq.n_kmers);
  }

  if( n_threads < freq.n_seq ) {

    pthread_mutex_init(&lock,NULL);
    
    kmers_t = (char***)malloc(sizeof(char**)*n_threads);
  
    for(i=0;i<n_threads;i++) {
      
      kmers_t[i] = gen_dna_kmers(kmer_length);
      
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
      free_kmers(kmers_t[i],kmer_length);
      pthread_mutex_destroy(&lock);
    }

    stride = freq.n_seq/n_threads;
    
    for(i=stride*n_threads;i<freq.n_seq%n_threads;i++) {
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
  
  free_kmers(kmers,kmer_length);
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

