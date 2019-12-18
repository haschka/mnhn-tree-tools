#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<unistd.h>

#include"dataset.h"
#include"cluster.h"

#if defined(_SCAN_SMITH_WATERMAN_GPU)
#ifdef __APPLE__
#include<OpenCL/OpenCL.h>
#else
#include<CL/opencl.h>
#endif
#endif

#include"dbscan.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {
  
  dataset ds;

  float epsilon;
  int minpts;

  FILE* fasta_f;
  
  split_set set_of_clusters;

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  opencl_stuff ocl;
#endif  

  if(argc < 4) {
    printf("Arguments are: \n"
	   " [file] Sequences in FASTA \n"
	   " Epsilon \n"
	   " minPoints \n"
	   " output-files prefix \n" );
    return(1);
  }

  sscanf(argv[2], "%f", &epsilon);
  sscanf(argv[3], "%i", &minpts);

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  ocl = opencl_initialization(ds);
  set_of_clusters = dbscan_SW_GPU(ds, epsilon, minpts, ocl);
#else
  set_of_clusters = dbscan_SW(ds, epsilon, minpts);
#endif
  
  printf("%u clusters obtained \n", set_of_clusters.n_clusters);

  if (set_of_clusters.n_clusters < 500) {
    create_cluster_files(argv[4], set_of_clusters, ds);
  }
  free_sequences_from_dataset(ds);
}
  
