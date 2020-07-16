#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<unistd.h>
#include"dataset.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

void print_arguments() {

  printf("Arguments are: \n"
	 "   [file] FASTA file to obtain sequence lengths from \n");
}

int main(int argc, char** argv) {

  dataset ds;

  int i;
  
  FILE* f = fopen(argv[1],"r");

  if (argc < 2) {
    print_arguments();
    return(1);
  }
 
  if ( NULL == f ) file_error(argv[1]);
  
  ds = dataset_from_fasta(f);
  
  for(i=0;i<ds.n_values;i++) {
    printf("%llu\n",(long long unsigned int)ds.sequence_lengths[i]);
  }
  
  free_sequences_from_dataset(ds);
  return(0);
}
