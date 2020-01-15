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

  int n_split_sets = argc - 4; 

  if (argc < 5) {
    printf("Arguments are: \n"
	   "  [fasta] file whose first sequence is the sequence to be found\n"
	   "  [fasta] database of sequences the split sets are refering to\n"
	   "  [string] color CSS compatible i.e. red, green, blue, #1122FF\n"
	   "  [split_sets...] cluster files corresponding to newick tree\n");
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

  printf("stroke:%s I ",argv[3]);
  for(i = 0;i<n_split_sets;i++) {
    s = read_split_set(argv[i+4]);
    for(j = 0;j<s.n_clusters; j++) {
      intersection = intersection_of_clusters(s.clusters[j], cl_search);
      if (intersection.n_members > 0) {
	printf("L%iC%iN%i ",i,j,s.clusters[j].n_members);
	free(intersection.members);
      }
    }
  }
  printf("\n");
  
  free_sequences_from_dataset(ds);
  free_sequences_from_dataset(ds_search);
  free(cl_search.members);
  return(0);
}
    
  
