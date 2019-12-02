#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<unistd.h>

#include"dataset.h"
#include"cluster.h"
#include"dbscan.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {
  
  dataset ds;

  size_t dimensions;
  float epsilon;
  int minpts;
  
  FILE* fasta_f;
  FILE* projection_f;
  
  split_set set_of_clusters;

  if(argc < 6) {
    printf("Arguments are: \n"
	   " [file] Sequences in FASTA \n"
	   " [file] Corresponding PCA from kmers \n"
	   " dimensions: dimensions in projections file"
	   " Epsilon \n"
	   " minPoints \n"
	   " output-files prefix \n" );
    return(1);
  }

  sscanf(argv[3], "%lu", &dimensions);
  sscanf(argv[4], "%f", &epsilon);
  sscanf(argv[5], "%i", &minpts);

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (fasta_f = fopen(argv[2], "r"))) file_error(argv[2]);
  projection_f = fopen(argv[2], "r");
  load_projections_from_file_into_dataset(projection_f,dimensions,&ds);
  fclose(projection_f);

  set_of_clusters = dbscan_L2(ds, epsilon, minpts);

  printf("%u clusters obtained \n", set_of_clusters.n_clusters);

  if (set_of_clusters.n_clusters < 500) {
    create_cluster_files(argv[6], set_of_clusters, ds);
  }
}
  
