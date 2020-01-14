#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include"dataset.h"
#include"cluster.h"


void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  int i, j;
  
  FILE* fasta_s;
  FILE* fasta_f;

  dataset ds_search;
  dataset ds;
  split_set s;

  cluster cl_search;
  cluster intersection;

  int n_split_sets = argc - 3; 

  if (argc < 4) {
    printf("Arguments are: \n"
	   "  [fasta] file whose first sequence is the sequence to be found\n"
	   "  [fasta] database of sequences the split sets are refering to\n"
	   "  [split_sets...] cluster files \n");
  }
  
  if ( NULL == (fasta_s = fopen(argv[1], "r"))) file_error(argv[1]);
  ds_search = dataset_from_fasta(fasta_s);
  fclose(fasta_s);

  if ( NULL == (fasta_f = fopen(argv[2], "r"))) file_error(argv[2]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  cl_search = cluster_from_sequence_in_dataset(ds, ds_search.sequences[0],
					       ds_search.sequence_lengths[0]);

  if (cl_search.n_members == 0) {
    printf("Sequence not in dataset \n");
    return(1);
  }

  printf("Sequence has internal ids: \n");

  for(i=0;i<cl_search.n_members;i++) {
    printf("%i ",cl_search.members[i]);
    if (i != 0 && i%8 == 0) printf("\n");
  }
  printf("\n");
  
  printf("Sequence is found in Clusters\n");
  
  for(i = 0;i<n_split_sets;i++) {
    s = read_split_set(argv[i+3]);
    for(j = 0; j<s.n_clusters; j++) {
      intersection = intersection_of_clusters(s.clusters[j], cl_search);
      if (intersection.n_members > 0) {
	printf("L%iC%i %i times\n",i,j,intersection.n_members);
	free(intersection.members);
      }
    }
  }
  
  free_sequences_from_dataset(ds);
  free_sequences_from_dataset(ds_search);
  free(cl_search.members);
  return(0);
}
    
  
