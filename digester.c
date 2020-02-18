#include<stdio.h>
#include<unistd.h>
#include<stddef.h>

#include"dataset.h"
#include"restriction_digest.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  dataset ds;
  FILE* fasta_f;

  dataset result;

  int n_threads;

  sscanf(argv[2],"%i",&n_threads);
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

#if defined(_digest_XbaI)  
  result = digest_XbaI(ds.sequences[0],ds.sequence_lengths[0],n_threads);
#elif defined(_digest_XmnI)
  result = digest_XmnI(ds.sequences[0],ds.sequence_lengths[0],n_threads);
#elif defined(_digest_HindIII)
  result = digest_HindIII(ds.sequences[0],ds.sequence_lengths[0],n_threads);
#endif
  dataset_to_fasta(stdout,result);
}
