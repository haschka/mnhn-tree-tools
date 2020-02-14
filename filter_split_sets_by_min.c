#include<stdio.h>
#include<stdlib.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"dataset.h"
#include"cluster.h"

int main(int argc, char** argv) {

  int i;
  int n_sets = argc-2;
  split_set *s = (split_set*)malloc(sizeof(split_set)*argc);

  split_set *sf;

  char path[256];
  
  if(argc < 3) {
    printf("Arguments are: "
	   " [int] minimum to consider a cluster \n"
	   " [path/prefix] Where to store the filtered splitsets\n "
	   " [file1 ... filen] split_sets obtainted from adaptive clustering \n"
           "                   calculations in order from lowest to hight\n");
  }

  sscanf(argv[1], "%i", &min_members);
  sscanf(argv[2], "%s", &path);
  
  for(i=3;i<argc;i++) {
    s[i-3] = read_split_set(argv[i]);
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

  root = generate_tree(n_sets, connections, s);
  print_tree(stdout, root);
}
