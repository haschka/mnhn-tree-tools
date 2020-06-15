#include<stdio.h>
#include<stdlib.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<unistd.h>

#include"dataset.h"
#include"cluster.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  int i;
  int clusters_offset;
  int n_sets = argc-2;
  split_set *s = (split_set*)malloc(sizeof(split_set)*argc);

  int lengths_from_outputfile;
  
  split_set *sf;
  
  tree_node* root;

  double* epsilons;
  double* densities;
  double* lengths;

  int min_members;
  int min_pts;

  int dimensions;

  FILE* cluster_out_f;
  
  cluster_connections** connections =
    (cluster_connections**)malloc(sizeof(cluster_connections*)*argc);

  if(argc < 2) {
    printf("Arguments are: \n"
	   " [int] filter for small clusters 0 = no filter \n"
	   "                                 n = minimum of n members in"
	   " cluster \n"
	   " [int] lengths from outputfile 0 = no file \n"
	   "                               1 = L1-density (min_pts/hypercube)\n"
	   "                               "
	   "2 = L2-density (min_pts/hypersphere)\n"
	   " [int] dimensions for density omitted if lengths = 0 \n"
	   " [outfile] clustering run outfile if lengths = (1 or 2) \n"
	   " [int] points for density calculations i.e. what min_pts you used\n"
	   " [file1 ... filen] split_sets obtainted from adaptive clustering \n"
           "                   calculations in order from lowest to hight\n");
    return(1);
  }

  sscanf(argv[1], "%i", &min_members);
  sscanf(argv[2], "%i", &lengths_from_outputfile);
  sscanf(argv[3], "%i", &dimensions);
  sscanf(argv[5], "%i", &min_pts);

  switch(lengths_from_outputfile) {
  case 0:
    clusters_offset = 3;
    break;
  default:
    clusters_offset = 6;
    break;
  }

  n_sets = argc - clusters_offset;
  
  for(i=clusters_offset;i<argc;i++) {
    s[i-clusters_offset] = read_split_set(argv[i]);
  }

  if(lengths_from_outputfile) {
    if ( NULL == (cluster_out_f = fopen(argv[4], "r"))) file_error(argv[4]);

    epsilons =
      get_epsilon_dist_from_adaptive_clustering_output(cluster_out_f,n_sets);
    fclose(cluster_out_f);
    densities = densities_from_epsilons(lengths_from_outputfile, dimensions,
					epsilons, n_sets, min_pts);
    free(epsilons);

    lengths = array_deltas(densities, n_sets);
    free(densities);
  }

  
  if(min_members) {
    sf = (split_set*)malloc(sizeof(split_set)*argc);
    for(i=0;i<n_sets;i++) {
      sf[i] = filtered_split_set_by_min_size(s[i],min_members);
      free_split_set_and_associated_clusters(s[i]);
      s[i] = sf[i];
    }
  }

  for(i=1;i<n_sets;i++) {
    connections[i-1] = generate_split_set_relation(s[i-1],s[i]);
  }

  switch(lengths_from_outputfile) {
  case 0:
    root = generate_tree(n_sets, connections, s, NULL);
    break;
  default:
    root = generate_tree(n_sets, connections, s, lengths);
    break;
  }

  for(i=0;i<n_sets;i++) {
    free_split_set_and_associated_clusters(s[i]);
  }
  
  free(s);
  print_tree(stdout, root);
}
