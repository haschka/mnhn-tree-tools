#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<unistd.h>
#include"binary_array.h"
#include"dataset.h"


void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  int i;
  
  dataset ds;

  FILE* fasta_f;
  FILE* mask_f;

  int i_buffer;

  char* binary_mask;

  if(argc < 3) {
    printf("Arguments are:\n");
    printf(" [file] fasta file to load sequences from\n");
    printf(" [file] mask - file containing\n"
	   "               \"1\"s for sequences to be inverted\n"
	   "               \"2\"s for sequences to be conserved\n");
    return(1);
  }
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (mask_f = fopen(argv[2], "r"))) file_error(argv[1]);

  binary_mask = alloc_and_set_zero_binary_array((size_t)ds.n_values);
  
  for(i = 0; i < ds.n_values; i++) {
    if (1 != fscanf(mask_f,"%i",&i_buffer)) {
      printf("Mask file currupt!\n");
      return(1);
    }
    if(i_buffer == 1) {
      set_value_in_binary_array_at_index(binary_mask,i);
    }
  }

  reverse_sequences(&ds, binary_mask);

  dataset_to_fasta(stdout, ds);

  free(binary_mask);
  free_sequences_from_dataset(ds);

  return(0);
  
}
      
    
  
