#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include"dataset.h"
#include"cluster.h"


void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  int i,j,k,l, count;
  
  char* found;
  int n_partitions;
  int n_layers;
  int partition_size;

  split_set set;
  
  FILE* fasta_f;

  int n_members ;
  int* members;

  float precision;

  dataset ds; 
  
  if (argc < 4) {
    printf(" Arguments are:\n"
	   " [FASTA] simulated dataset of n partitions\n"
	   " [int] n number of partitions in the dataset\n"
	   " [float] precision \n"
	   " [split_sets] partition files from a clustering run\n"
	   );
    return(1);
  }
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  sscanf(argv[2], "%i", &n_partitions);
  sscanf(argv[3], "%f", &precision);

  if(precision > 1. || precision < 0. ) {
    printf("Error ! requested cluster presion not in range [0-1]");
  }
  
  precision = 1. - precision;
  
  if(ds.n_values%n_partitions != 0 ) {
    printf("Error ! Number of sequences is not divisiable by the number of \n"
	   "partitions \n");
  }
  
  partition_size = ds.n_values/n_partitions;
  
  n_layers = argc-4;

  found = (char*)malloc(sizeof(char)*n_partitions);
  
  memset(found,0, sizeof(char)*n_partitions);

  for(i=0;i<n_layers;i++){
    set = read_split_set(argv[i+4]);

    for(j=0;j<n_partitions;j++) {
      for(k=0;k<set.n_clusters;k++) {
	count = 0;
	n_members = set.clusters[k].n_members;
	members = set.clusters[k].members;
	if((float)n_members <
	   (float)partition_size+(float)partition_size*precision &&
	   (float)n_members >
	   (float)partition_size-(float)partition_size*precision) {
	  for(l=0;l<n_members;l++) {
	    if(j*partition_size <= members[l] &&
	       members[l] < (j+1)*partition_size) {
	      count++;
	    }
	  }
	  if ((float)count >
	      (float)partition_size-(float)partition_size*precision) {
	    found[j] = 1;
	  }
	
	}
      }
    }
  }

  count = 0;
  for(i=0;i<n_partitions;i++) {
    if(found[i] == 1) { count++; }
  }
  printf("Found %i out of %i partitions\n",count, n_partitions);
  return(0);
}
  
  
    
  
