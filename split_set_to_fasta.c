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

  if(argc<4) {
    printf("Arguments are: \n"
	   "   1. [file] fasta file of complete dataset\n"
	   "   2. [file] file containing the cluster set\n"
	   "   3. path prefix where clusters should be stored to\n");
    return(1);
  }
  
  f=fopen(argv[1],"r");
  s=read_split_set(argv[2]);
  ds=dataset_from_fasta(f);
  
  create_cluster_files(argv[3],s,ds);
  return(0);
}

  
