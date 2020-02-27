#include<stdio.h>
#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"dataset.h"
#include"cluster.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  FILE* fasta_f; 
  FILE* projection_f;
  
  split_set s; 

  dataset ds;

  int dimensions;

  if(argc<4) {
    printf("Arguments are: \n"
	   "   1. [file] fasta file of complete dataset\n"
	   "   2. [file] pca projections for entire fasta file\n"
	   "   3. (int) number of principal components in projections file\n"
	   "   4. [file] split set containing the cluster set\n"
	   "   5. path prefix where clusters should be stored to\n");
    return(1);
  }

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  sscanf(argv[3],"%i",&dimensions);

  
  if ( NULL == (projection_f = fopen(argv[2], "r"))) file_error(argv[2]);
  load_projections_from_file_into_dataset(projection_f,dimensions,&ds);
  fclose(projection_f);

  
  s=read_split_set(argv[4]);
  
  create_cluster_files_with_values(argv[5],s,ds);
  return(0);
}

  
