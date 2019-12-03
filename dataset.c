#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<unistd.h>

#include"dataset.h"

void char_sequence_to_binary(char* c_seq, char* b_seq, size_t seq_len) {

  size_t i;
  int j;
  char bin;
  for(i = 0; i< seq_len; i++) {

    switch(c_seq[i]) {
    case 'A':
      bin = 0;
      break;
    case 'C':
      bin = 1;
      break;
    case 'G':
      bin = 2;
      break;
    case 'T':
      bin = 3;
      break;
    default:
      printf("unhandled character for binary conversion %c", c_seq[i]);
      _exit(1);
    }

    j=i%4;

    b_seq[i/4] |= (bin << 2*j);
  }
}

char get_char_from_binary_sequence_at_index(char* b_seq, size_t idx) {

  size_t b_idx = idx/4;
  size_t b_off = idx%4;

  char test_a = (1 << (b_off*2));
  char test_b = (1 << (b_off*2)+1);

  char res = 0;

  char ret_val;
  
  if(test_a == (b_seq[b_idx] & test_a)) res += 1;
  if(test_b == (b_seq[b_idx] & test_b)) res += 2;

  switch (res) {
    case 0:
      ret_val = 'A';
      break;
    case 1:
      ret_val = 'C';
      break;
    case 2:
      ret_val = 'G';
      break;
    case 3:
      ret_val = 'T';
      break;
    }
  return(ret_val);
}
  

void add_binary_sequences_to_dataset(dataset* ds) {

  size_t i;
  char** c_seq = ds->sequences;
  size_t* lengths = ds->sequence_lengths;
  char** b_seq = malloc(sizeof(ds->n_values)*sizeof(char*));
  
  for(i = 0 ; i< ds->n_values; i++ ) {
    b_seq[i] = (char*)malloc(sizeof(char)*(lengths[i]/4+1));
    char_sequence_to_binary(c_seq[i], b_seq[i], lengths[i]);
  }
}

dataset dataset_from_fasta(FILE* in) {

  char* line = NULL;
  size_t line_size = 0;

  int sequences = 0;

  dataset ds;

  char linebuffer[2000];
  size_t linebuffer_length;
  size_t sequence_length;
  int i;
  
  rewind(in);

  while ( -1 != getline(&line, &line_size, in) ) {
    if( line[0] == '>' ) sequences++;
  }

  rewind(in);

  ds.sequence_lengths = (size_t*)malloc(sizeof(size_t*)*sequences);

  ds.sequences = (char**)malloc(sizeof(char*)*sequences);
  ds.n_values = sequences;

  sequences = -1;

  ds.max_sequence_length = 0;
  
  while ( -1 != getline(&line, &line_size, in) ) {
    
    if ( line[0] == '>') {
      /* new sequence in file */
      sequences++;

      /* redress memory usage of previous sequence */
      if(sequences > 0) {
	ds.sequences[sequences-1] =
	  (char*)realloc(ds.sequences[sequences-1],
			 sizeof(char)*(ds.sequence_lengths[sequences-1]+1)); 
      }
			 						   
      ds.sequences[sequences] = (char*)malloc(sizeof(char)*1001);
      ds.sequences[sequences][0] = 0;
      ds.sequence_lengths[sequences] = 0;

    } else {
      /* continuation of current sequence */

      sscanf(line,"%s",linebuffer);

      linebuffer_length = strlen(linebuffer);
      ds.sequence_lengths[sequences] += linebuffer_length;
      if(ds.sequence_lengths[sequences] > 1000) {
	ds.sequences[sequences] =
	  (char*)realloc(ds.sequences[sequences],
			 sizeof(char)*(ds.sequence_lengths[sequences]+1));
      }
      for(i=0; i < linebuffer_length; i++) {
	linebuffer[i] = toupper(linebuffer[i]);
      }
      strcat(ds.sequences[sequences],linebuffer);
    }
    
  }
  for(i=0;i<ds.n_values;i++) {
    if(ds.max_sequence_length < ds.sequence_lengths[i]) {
      ds.max_sequence_length = ds.sequence_lengths[i];
    }
  }
  free(line);
  return(ds);
}

data_shape shape_from_kmer_file(int infile) {

  data_shape s;

  unsigned char current_character;

  off_t size;

  int i;

  s.n_features = 0;
  s.n_samples = 0;

  size = lseek(infile, 0, SEEK_END);
  lseek(infile, 0, SEEK_SET);

  char* f_buffer = (char*)malloc(sizeof(char)*size);

  if (size != read(infile, f_buffer, size)) {
    printf("Could not read file in order to read obtain data shape\n");
    _exit(1);
  }

  i = 0;
  while( f_buffer[i] != '\n') {
    if ( f_buffer[i] == '\t') s.n_features++;
    i++;
  }

  while( i < size ) {
    if ( f_buffer[i] == '\n' ) s.n_samples++;
    i++;
  }

  return(s);
}


dataset load_kmer_from_file_into_dataset(FILE* in_file, data_shape shape) {

  dataset ds;
  int i, j;

  char buffer[1024];

  ds.n_dimensions = shape.n_features;

  ds.n_values = shape.n_samples;

  ds.values = (float**)malloc(sizeof(float*)*ds.n_dimensions);

  for(i = 0 ; i < ds.n_dimensions; i++) {

#if defined(__AVX__) || defined(__SSE__)
    if(0 != posix_memalign((void**)&ds.values[i],
#if defined(__AVX__)
			   32,
#elif defined(__SSE__)
			   16,
#endif
			   sizeof(float)*ds.n_values)) {
      printf("Error allocating aligned memory for kmers \n");
      _exit(1);
    }
#else
    ds.values[i] =
      (float*)malloc(sizeof(float)*ds.n_values);
#endif
  }

  rewind(in_file);

  for(i=0;i<shape.n_samples;i++) {
    fscanf(in_file,"%s", buffer);
    for(j=0;j<shape.n_features;j++) {
      fscanf(in_file,"%f", ds.values[j]+i);
    }
  }
  return(ds);
}


void load_projections_from_file_into_dataset(FILE* projections,
					     size_t dimensions,
					     dataset* ds) {

  size_t i,j;

  size_t counter;

  ds->n_dimensions = dimensions;
  ds->values = (float**)malloc(sizeof(float*)*ds->n_dimensions);

  for(i = 0 ; i < ds->n_dimensions; i++) {

#if defined(__AVX__) || defined(__SSE__)
    if(0 != posix_memalign((void**)&ds->values[i],
#if defined(__AVX__)
			   32,
#elif defined(__SSE__)
			   16,
#endif
			   sizeof(float)*(ds->n_values))) {
      printf("Could not allocate memory to load projections \n");
      _exit(1);
    }
#else
    if (NULL == (ds->values[i] = (float*)malloc(sizeof(float)*ds.n_values))) {
      printf("Could not allocate memory to load projections \n");
      _exit(1);
    }    
#endif
  }

  for(i=0;i<ds->n_values;i++) {
    counter = 0;
    for(j=0;j<dimensions;j++) {
      counter += fscanf(projections, "%f", ds->values[j]+i);
    }
    if(counter != dimensions) {
      printf("Error reading projections file \n");
      _exit(1);
    }
  }
}



void free_values_from_dataset(dataset ds) {
  int i;
  for(i=0;i<ds.n_dimensions;i++) {
    free(ds.values[i]);
  }
  free(ds.values);
}

void free_sequences_from_dataset(dataset ds) {
  int i;
  for(i=0;i<ds.n_values;i++) {
    free(ds.sequences[i]);
  }
  free(ds.sequences);
}

void free_dataset(dataset ds) {
  free_values_from_dataset(ds);
  free_sequences_from_dataset(ds);
}
