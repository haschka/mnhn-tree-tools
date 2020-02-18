#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include"dataset.h"
#include"filter.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc,char** argv) {

  FILE* fasta_f;

  FILE* target_f;

  dataset ds, target_set;
  dataset new;

  size_t size,plus,minus;
  unsigned int max_distance;
  int n_threads;

  if (argc < 7) {
    printf("Arguments are: \n"
	   " [file] Fasta dataset to filter \n"
	   " [file] Fasta target sequence for Smith Waterman filter \n"
	   " (int) size filter size \n"
	   " (int) allow sequences in the range [size size+this_value] \n"
	   " (int) allow sequences in the tange [size-this_value size] \n"
	   " (int) maximum allowed Smith Waterman distance from target \n"
	   " (int) number of threads to use for this computation \n");
    return(1);
  }

  sscanf(argv[3],"%lu",&size);
  sscanf(argv[4],"%lu",&plus);
  sscanf(argv[5],"%lu",&minus);
  sscanf(argv[6],"%u",&max_distance);
  sscanf(argv[7],"%i",&n_threads);

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (target_f = fopen(argv[2], "r"))) file_error(argv[2]);
  target_set = dataset_from_fasta(target_f);
  fclose(target_f);

  new = filter_ds_by_size(ds,size,plus,minus);
  free_sequences_from_dataset(ds);
  ds = new;

  new = filter_by_SW_distance(ds,
			      target_set.sequences[0],
			      target_set.sequence_lengths[0],
			      max_distance,
			      n_threads);
  free_sequences_from_dataset(ds);
  ds = new;

  dataset_to_fasta(stdout,ds);
  free_sequences_from_dataset(ds);
  free_sequences_from_dataset(target_set);
  return(0);
}
  
				 
  
  
  
  

