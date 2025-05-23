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
	 "   [file] FASTA file to obtain unique sequences from \n");
}


int main(int argc, char** argv) {

  dataset ds;

  FILE* input = fopen(argv[1],"r");

  unique_sequences us;

  if (argc < 2) {
    print_arguments();
    return(1);
  }
  
  if (input == NULL) file_error(argv[1]);
  
  ds = dataset_from_fasta(input);

  us = get_sequence_multiplicities(ds);
  sort_unique_sequences(us);
  write_unique_sequences(stdout, ds, us);

  fclose(input);
  
}
  
