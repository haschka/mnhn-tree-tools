#include<stdio.h>
#include"dataset.h"  

int main(int argc, char** argv) {

  FILE* f;

  dataset ds;
  consens cs;
  
  if(argc < 2) {
    printf("Argument: [file] Fasta file to built consens from \n");
    return(0);
  }

  f = fopen(argv[1],"r");

  ds = dataset_from_fasta(f);          
  cs = obtain_consens_from_dataset(ds);
  
  print_consensus_statistics(stdout,cs);
}
  
  
  
