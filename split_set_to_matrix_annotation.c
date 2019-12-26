#include<stdio.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"dataset.h"
#include"cluster.h"

int main(int argc, char** argv) {

  FILE* f; 

  dataset ds; 

  if(argc<2) {
    printf("Arguments are: \n"
	   "   1. [file] fasta file of complete dataset\n");
    return(1);
  }
  
  f=fopen(argv[1],"r");
  ds=dataset_from_fasta(f);

  print_cluster_matrix_view_annotation(stdout,ds);
  return(0);
}

  
