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

  int i;
  
  FILE* fasta_f;

  dataset ds;
  dataset ds_target;

  split_set target;

  split_set *s;

  int clusters_offset = 3;
  int n_sets = argc - clusters_offset;

  cluster_connections** connections =
    (cluster_connections**)malloc(sizeof(cluster_connections*)*argc);

  tree_node* root;
  double* purness;
  double* cluster_factors;
  
  if (argc < 4) {
    printf("Arguments are: \n"
	   "  [fasta] database of sequences the split sets are refering to \n"
	   "  [split-set] target split set for purness evaluation \n"
	   "  [split_sets...] cluster files corresponding to a tree \n");
    return(1);
  }

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  s = (split_set*)malloc(sizeof(split_set)*n_sets);

  for(i=clusters_offset;i<argc;i++) {
    s[i-clusters_offset] = read_split_set(argv[i]);
  }

  target = read_split_set(argv[2]);
  
  for(i=1;i<n_sets;i++) {
    connections[i-1] = generate_split_set_relation(s[i-1],s[i]);
  }

  root = generate_tree(n_sets, connections, s, NULL);

  purness = pureness_from_tree(n_sets, root, s, target);
  cluster_factors = clusters_in_layer_vs_target_clusters(n_sets, s, target);
  
  for(i=0;i<n_sets;i++) {
    printf("%lf\t%lf\n",purness[i],cluster_factors[i]);
  }
}
