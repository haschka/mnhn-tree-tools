#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"dataset.h"
#include"binary_array.h"

/*! \brief a structure for unique sequences
 */
typedef struct {
  int* u_seq;            /*!< array with indexes to
                          *    unique sequences in a dataset */
  int* multiplicities;   /*!< array with the number of duplicates */
  int n_seq;             /*!< number of unique sequences in dataset */
} unique_sequences;

typedef struct {
  int buffer[2];
} sorthelper;

int comp(const void *a, const void *b) {

  sorthelper* sa = (sorthelper*)a;
  sorthelper* sb = (sorthelper*)b;

  return((sa->buffer[1]-sb->buffer[1]));
}

void sort_unique_sequences(unique_sequences us) {

  int i;
  sorthelper* helper = (sorthelper*)malloc(sizeof(sorthelper)*us.n_seq);

  for(i=0;i<us.n_seq;i++) {
    helper[i].buffer[0] = us.u_seq[i];
    helper[i].buffer[1] = us.multiplicities[i];
  }

  qsort(helper, us.n_seq, sizeof(sorthelper), comp);

  for(i=0;i<us.n_seq;i++) {
    us.u_seq[i] = helper[i].buffer[0];        
    us.multiplicities[i] = helper[i].buffer[1];
  }
}

unique_sequences get_sequence_multiplicities(dataset ds) {

  unique_sequences us;
  
  char* visited_sites = alloc_and_set_zero_binary_array(ds.n_values);

  int i,j, u_seq_counter;

  char* current_sequence;

  us.u_seq = (int*)malloc(sizeof(int)*ds.n_values);
  us.multiplicities = (int*)malloc(sizeof(int)*ds.n_values);

  for(i=0; i<ds.n_values;i++) {
    us.multiplicities[i] = 1;
  }
  
  u_seq_counter=0;
  for(i=0; i<ds.n_values;i++) {
    if (!get_value_in_binary_array_at_index(visited_sites,i)) {
      us.u_seq[u_seq_counter] = i;
      set_value_in_binary_array_at_index(visited_sites,i);
      for(j=0; j<ds.n_values;j++) {
	if(!get_value_in_binary_array_at_index(visited_sites,j)) {
	  if (!strcmp(ds.sequences[i],ds.sequences[j]) &&
	      (strlen(ds.sequences[i]) == strlen(ds.sequences[j]))) {
	    set_value_in_binary_array_at_index(visited_sites,j);
	    us.multiplicities[u_seq_counter]++;
	  }
	}
      }
      u_seq_counter++;
    }
  }

  us.u_seq = (int*)realloc(us.u_seq, sizeof(int)*u_seq_counter);
  us.multiplicities =
    (int*)realloc(us.multiplicities, sizeof(int)*u_seq_counter);
  us.n_seq = u_seq_counter;
  return(us);
}

void write_unique_sequences(FILE* outfile, dataset ds, unique_sequences us) {

  int i,k;
  for(i=0; i<us.n_seq; i++) {
    fprintf(outfile,">sequence_%i_mutiplicity_%i\n",
	    us.u_seq[i],
	    us.multiplicities[i]);
    for( k = 0;
	 k < (strlen(ds.sequences[us.u_seq[i]]))-1;
	 k++) {
      if( k != 0 && k%50 == 0 ) {
	fprintf(outfile, "\n");
	fputc(ds.sequences[us.u_seq[i]][k],outfile);
      } else {
	fputc(ds.sequences[us.u_seq[i]][k],outfile);
      }
    }
    if( (strlen(ds.sequences[us.u_seq[i]])-1) != 0 &&
	(strlen(ds.sequences[us.u_seq[i]])-1) %50 == 0 ) {
      fprintf(outfile, "\n");
      fputc(ds.sequences[us.u_seq[i]][k], outfile);
    } else {
      fputc(ds.sequences[us.u_seq[i]][k], outfile);
    }
    fprintf(outfile,"\n");
  }
}	      
  
int main(int argc, char** argv) {

  dataset ds;

  FILE* input = fopen(argv[1],"r");

  unique_sequences us;
  
  ds = dataset_from_fasta(input);

  us = get_sequence_multiplicities(ds);
  sort_unique_sequences(us);
  write_unique_sequences(stdout, ds, us);

  fclose(input);
  
}
  
