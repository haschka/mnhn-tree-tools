#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<unistd.h>

#if defined(_SCAN_SMITH_WATERMAN_GPU)
#ifdef __APPLE__
#include<OpenCL/OpenCL.h>
#else
#include<CL/opencl.h>
#endif
#endif

#include"dataset.h"
#include"cluster.h"
#include"dbscan.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  int i,j,k;
  
  dataset ds;

  float epsilon;
  int minpts;

  FILE* fasta_f;

  float epsilon_start;
  float epsilon_inc;

  int initial_counter, count, eps_count;

  char split_files[255];
  char split_files_prefix[255];
  char buffer[20];
  
  split_set* set_of_split_sets;
  cluster_connections** connections;
  cluster_connections* current_connection;
  
  cluster not_covered;

  float coverage;
  
  split_set new_split_set;

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  opencl_stuff ocl;
#endif
  
  sscanf(argv[2], "%f", &epsilon_start);
  sscanf(argv[3], "%f", &epsilon_inc);
  sscanf(argv[3], "%i", &minpts);
  sscanf(argv[4], "%s", split_files_prefix);
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  ocl = opencl_initialization(ds);
#endif
  
  set_of_split_sets = (split_set*)malloc(sizeof(split_set));

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  set_of_split_sets[0] = dbscan_SW_GPU(ds, epsilon_start, minpts,ocl);
#else
  set_of_split_sets[0] = dbscan_SW(ds, epsilon_start, minpts);
#endif
 

  initial_counter = 0;
  while( set_of_split_sets[0].n_clusters == 1 ) {
    if (initial_counter == 20) {
      printf("Error did not find a sufficient" 
	     " starting position in 20 tries \n");
#if defined(_SCAN_SMITH_WATERMAN_GPU)
      set_of_split_sets[0] = dbscan_SW_GPU(ds, epsilon_start/2, minpts,ocl);
#else
      set_of_split_sets[0] = dbscan_SW(ds, epsilon_start/2, minpts);
#endif
      initial_counter++;
    }  
  }

  printf("Sarting with %i clusters\n", set_of_split_sets[0].n_clusters);
  
  not_covered = data_not_in_clusters(set_of_split_sets[0], ds);

  printf("Coverage at initial point: %f",
	 (float)1.f-(float)not_covered.n_members/(float)ds.n_values);

  if(set_of_split_sets[0].n_clusters == 1) {
    printf("All clusters fusioned in one step, decrease epsilon increment\n");
    return(1);
  }

  count = 0;
  eps_count = 1;
  do{
#if defined(_SCAN_SMITH_WATERMAN_GPU)
    new_split_set =
      dbscan_SW_GPU(ds,epsilon_start+eps_count*epsilon_inc, minpts, ocl);
#else
    new_split_set = dbscan_SW(ds, epsilon_start+eps_count*epsilon_inc, minpts);
#endif
    if (new_split_set.n_clusters < set_of_split_sets[count].n_clusters) {
      
      not_covered = data_not_in_clusters(new_split_set, ds);
      
      printf("Coverage at layer %i: %f", count,
	     1.f-(float)not_covered.n_members/(float)ds.n_values);
	          
      current_connection = generate_split_set_relation(set_of_split_sets[0],
						       new_split_set);
      count++;
      set_of_split_sets = (split_set*)realloc(set_of_split_sets,
					      sizeof(split_set)*(count+1));
      connections =
	(cluster_connections**)realloc(connections,
				       sizeof(cluster_connections*)*
				       count);
      
      set_of_split_sets[count] = new_split_set;
      connections[count-1] = current_connection;
    }
    eps_count++;
  }while(new_split_set.n_clusters != 1);

  printf("Connections found: \n");

  for(i=0;i<count;i++) {
    
    sprintf(buffer,"%04d",i);
    memcpy(split_files,split_files_prefix,
	   strlen(split_files_prefix)+1);
    strcat(split_files,buffer);

    store_split_set(split_files, set_of_split_sets[i]);
  }

  for(i=count;i>0;i--) {
    printf("Layer %i:\n", i);
    for(j=0;set_of_split_sets[i].n_clusters;j++) {
      printf("  Cluster %i connected to: \n  ",j);
      for(k=0;connections[i-1][j].n_connections;k++) {
	printf("%i ", connections[i-1][j].connections[k]);
      }
    }
  }
}
  
  
  

    
    
    
  
    
