#include<stdio.h>
#include<stdlib.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"dataset.h"
#include"cluster.h"

int main(int argc, char** argv) {

  int i;
  int n_sets = argc-1;
  split_set *s = (split_set*)malloc(sizeof(split_set)*argc);

  tree_node* root;
  
  cluster_connections** connections =
    (cluster_connections**)malloc(sizeof(cluster_connections*)*argc);
  
  for(i=1;i<argc;i++) {
    s[i-1] = read_split_set(argv[i]);
  }

  for(i=1;i<n_sets;i++) {
    connections[i-1] = generate_split_set_relation(s[i-1],s[i]);
  }

  root = generate_tree(n_sets, connections, s);
  print_tree(stdout, root);
}
