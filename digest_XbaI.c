#include<stdio.h>
#include<unistd.h>
#include<stddef.h>

#include"dataset.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

dataset digest_XbaI(char * sequence, size_t sequence_length);

int main(int argc, char** argv) {

  dataset ds;
  FILE* fasta_f;

  dataset result;

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  result = digest_XbaI(ds.sequences[0],ds.sequence_lengths[0]);
  
  dataset_to_fasta(stdout,result);
}
