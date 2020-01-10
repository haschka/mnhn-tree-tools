#include<stdio.h>

#include"dataset.h"
#include"kmers.h"

int main(int argc, char** argv) {

  size_t length;
  dataset ds;
  kmer_frequencies freq;
  int i;
  int n_threads;
  
  FILE* in = fopen(argv[1],"r");

  if(argc < 3) {
    printf("Arguments are: \n"
	   "  [file] fasta file \n"
	   "  [int] kmer length \n"
	   "  [int] number of threads \n");
    return(1);
  }
  
  sscanf(argv[2],"%lu", &length);
  sscanf(argv[3],"%d", &n_threads);
  ds = dataset_from_fasta(in);
  fclose(in);
  freq = frequencies_from_dataset(ds, length, n_threads);
  write_kmer_base(stdout, freq);

  free_kmer_frequencies(freq);
  free_sequences_from_dataset(ds);
    
}
