#include<stdio.h>
#include"dataset.h"  

int main(int argc, char** argv) {

  FILE* f = fopen(argv[1],"r");

  dataset ds = dataset_from_fasta(f);
  consens cs = obtain_consens_from_dataset(ds);

  print_consensus_statistics(stdout,cs);
}
  
  
  
