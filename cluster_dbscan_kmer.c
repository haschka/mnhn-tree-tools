#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include"dataset.h"
#include"cluster.h"
#include"dbscan.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {
  
  dataset ds;

  float epsilon;
  int minpts;

  data_shape shape;
  
  FILE* fasta_f;
  int kmer_fd;
  FILE* kmer_f;
  
  split_set set_of_clusters;

  if(argc < 5) {
    printf("Arguments are: \n"
	   " [file] Sequences in FASTA \n"
	   " [file] Corresponding kmers \n"
	   " Epsilon \n"
	   " minPoints \n"
	   " output-files prefix \n" );
    return(1);
  }

  sscanf(argv[3], "%f", &epsilon);
  sscanf(argv[4], "%i", &minpts);

  if ( -1 == (kmer_fd = open(argv[2], O_RDONLY))) file_error(argv[2]);
  
  shape = shape_from_kmer_file(kmer_fd);
  kmer_f = fdopen(kmer_fd,"r");
  ds = load_kmer_from_file_into_dataset(kmer_f, shape);
  fclose(kmer_f);
#if defined(_CLUSTER_KMER_L1)
  set_of_clusters = dbscan_L1(ds, epsilon, minpts);
#elif defined(_CLUSTER_KMER_L2)
  set_of_clusters = dbscan_L2(ds, epsilon, minpts);
#endif

  free_values_from_dataset(ds);
  
  printf("%u clusters obtained \n", set_of_clusters.n_clusters);

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);
  if (set_of_clusters.n_clusters < 500) {
    create_cluster_files(argv[5], set_of_clusters, ds);
  }
  free_sequences_from_dataset(ds);
}
  
