#include<stdio.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"dataset.h"
#include"cluster.h"

int main(int argc, char** argv) {

  FILE* f; 

  split_set s; 

  dataset ds; 

  if(argc<3) {
    printf("Arguments are: \n"
	   "   1. [file] fasta file of complete dataset\n"
	   "   2. [file] file containing the cluster set\n");
    return(1);
  }
  
  f=fopen(argv[1],"r");
  s=read_split_set(argv[2]);
  ds=dataset_from_fasta(f);
  
  print_cluster_matrix_view(stdout, s, ds);
  return(0);
}

  
