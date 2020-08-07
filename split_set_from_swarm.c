#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#include"dataset.h"
#include"cluster.h"
#include"binary_array.h"

void file_error(char* path) {
  printf("failed to open file %s\n", path);
  _exit(1);
}

void print_arguments() {

  printf("Arguments are: \n");
  printf(" [FILE-in]  Fasta file to cluster according to swarm file \n");
  printf(" [FILE-in]  Swarm file \n");
  printf(" [FILE-out] split-set file to be created \n");
  _exit(1);
}

int main(int argc, char** argv) {

  int i;
  
  dataset ds;

  FILE* fasta_f;
  FILE* swarm_f;

  split_set s;

  char buffer[1000];

  char * visited_sites_b;

  cluster * current_cluster;
  size_t length;

  int sequence_no, multiplicity;
  
  if(argc < 3) print_arguments();

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (swarm_f = fopen(argv[2], "r"))) file_error(argv[2]);
  
  visited_sites_b = alloc_and_set_zero_binary_array(ds.n_values);

  s.n_clusters = 0;
  s.clusters = (cluster*)malloc(sizeof(cluster));
  s.clusters[0].members = (int*)malloc(sizeof(int));
  s.clusters[0].n_members = 0;
  
  while ( 1 == fscanf(swarm_f,"%s",buffer) ) {
    
    sscanf(buffer,"sequence_%i_mutiplicity_%i",&sequence_no, &multiplicity);
    s.clusters[s.n_clusters].members[s.clusters[s.n_clusters].n_members] =
      sequence_no;
    set_value_in_binary_array_at_index(visited_sites_b,sequence_no);
    s.clusters[s.n_clusters].n_members++;

    s.clusters[s.n_clusters].members = 
      (int*)realloc(s.clusters[s.n_clusters].members,
		    sizeof(int)*(s.clusters[s.n_clusters].n_members+1));
    
    if(multiplicity != 1) {

      s.clusters[s.n_clusters].members = 
	(int*)realloc(s.clusters[s.n_clusters].members,
		      sizeof(int)*(s.clusters[s.n_clusters].n_members
				   +multiplicity));

      for(i=0;i<ds.n_values;i++) {
	if(!get_value_in_binary_array_at_index(visited_sites_b,i)) {
	  if(!strcmp(ds.sequences[sequence_no],ds.sequences[i]) &&
	     (strlen(ds.sequences[sequence_no]) == strlen(ds.sequences[i]))) {
	    set_value_in_binary_array_at_index(visited_sites_b, i);

	    s.clusters[s.n_clusters].
	      members[s.clusters[s.n_clusters].n_members] = i;
	    
	    s.clusters[s.n_clusters].n_members++;
	  }
	}
      }
    }
    
    if(fgetc(swarm_f) == '\n') {
      s.n_clusters++;
      s.clusters = (cluster*)realloc(s.clusters,
				     sizeof(cluster)*(s.n_clusters+1));
      s.clusters[s.n_clusters].members = (int*)malloc(sizeof(int));
      s.clusters[s.n_clusters].n_members = 0;
    }
  }
  
  store_split_set(argv[3],s);
  free_sequences_from_dataset(ds);
  free_split_set_and_associated_clusters(s);
  return(0);
}
				     
	      
	   
		       
	    
  
