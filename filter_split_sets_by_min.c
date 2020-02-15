#include<stdio.h>
#include<stdlib.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<string.h>

#include"dataset.h"
#include"cluster.h"

int main(int argc, char** argv) {

  int i;
  int n_sets = argc-3;
  int min_members;
  split_set *s = (split_set*)malloc(sizeof(split_set)*argc);

  split_set *sf;

  char path[256];
  char number_in_file[10];

  char filename[256];

  FILE* f;
  
  filename[0] = 0;
  
  if(argc < 3) {
    printf("Arguments are: "
	   " [int] minimum to consider a cluster \n"
	   " [path/prefix] Where to store the filtered splitsets\n "
	   " [file1 ... filen] split_sets obtainted from adaptive clustering \n"
           "                   calculations in order from lowest to hight\n");
  }

  sscanf(argv[1], "%i", &min_members);
  sscanf(argv[2], "%s", path);
  
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

  for(i=0;i<n_sets;i++) {
    strcat(filename,path);
    sprintf(number_in_file,"%04d",i);
    strcat(filename,number_in_file);
    store_split_set(filename,s[i]);
    filename[0] = 0;
    free_split_set_and_associated_clusters(s[i]);
  }

}
