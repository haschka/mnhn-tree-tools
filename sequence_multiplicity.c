#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"dataset.h"
  
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
  
