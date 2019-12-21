#include<stdio.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"dataset.h"
#include"cluster.h"

int main(int argc, char** argv) {

  int i,j;
  
  FILE* f; 

  split_set s_a,s_b; 

  dataset ds; 

  cluster_connections* connections;
  
  if(argc<3) {
    printf("Arguments are: \n"
	   "   1. [file] lower_layer_split_set\n"
	   "   2. [file] upper_layer_split_set\n");
    return(1);
    }
  
  s_a=read_split_set(argv[1]);
  s_b=read_split_set(argv[2]);

  connections = generate_split_set_relation(s_a,s_b);

  for(i=0;i<s_b.n_clusters;i++) {
    printf("Cluster %i connected to \n", i);
    for(j=0;j<connections[i].n_connections;j++) {
      printf("%i ", connections[i].connections[j]);
    }
    printf("\n");
  }
  
  return(0);
}

  
