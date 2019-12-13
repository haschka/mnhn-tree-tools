#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<unistd.h>

#include"dataset.h"
#include"cluster.h"
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

  float epsilon_start;
  float epsilon_inc;
  int minpts;

  int count;
  
  split_set current_set_clusters;
  split_set old_set_clusters;

  cluster not_covered;
  
  sscanf(argv[2], "%f", &epsilon_start);
  sscanf(argv[3], "%f", &epsilon_inc);
  sscanf(argv[3], "%i", &minpts);
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  old_set_clusters = dbscan_SW(ds, epsilon_start, minpts);

  count = 0;
  while( old_set_clusters.n_clusters == 1 ) {
    if (counter == 20) {
      printf("Error did not find a sufficient" 
	     " starting position in 20 tries \n")'
    old_set_clusters = dbscan_SW(ds, epsilon_start/2, minpts);
    count++;
    }  
  }

  prinf("Sarting with %i clusters\n", old_set_clusters.n_clusters);
  
  not_covered = data_not_in_clusters(old_set_clusters, ds);

  printf("Coverage at initial point: %f",
	 (float)1.f-(float)not_covered/(float)ds.n_values);

  count = 1;
  

  if(new_set_clusters.n_clusters == 1) {
    printf("All clusters fusioned in one step, decrease epsilon increment\n");
    return(1);
  }

  do{
    new_set_clusters = dbscan_SW(ds, epsilon_start+count*epsilon_inc, minpts);
    if (new_set_clusters.n_clusters < old_set_clusters.n_clusters) {
      connections[count-1] = generate_split_set_relation(old_set_clusters,
							 new_set_clusters);
    }
	
  }while(new_set_clusters.n_cluster != 1)

    
    
    
  
    
