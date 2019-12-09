#include<stdio.h>
#include<string.h>
#include"dataset.h"
#include"cluster.h"

void create_cluster_files_with_values(char* prefix, split_set s, dataset ds) {

  int i,j,k;

  FILE* cur_file; 

  char file_name[256];

  char i_buffer[4];
  
  for ( i = 0; i< s.n_clusters; i++) {

    file_name[0] = 0;
    strcat(file_name, prefix);
    strcat(file_name, "-");
    sprintf(i_buffer,"%03d", i);
    strcat(file_name, i_buffer);

    cur_file = fopen(file_name, "w");

    for (j = 0; j< s.clusters[i].n_members; j++) {
      for ( k = 0; k < ds.n_dimensions - 1; k++) {
	fprintf(cur_file, "%f\t", ds.values[k][s.clusters[i].members[j]]);
      }
      fprintf(cur_file, "%f\n",
	      ds.values[ds.n_dimensions-1][s.clusters[i].members[j]]);
    }

    fclose(cur_file);
  }
}

void create_cluster_files(char* prefix, split_set s, dataset ds) {

  int i,j,k;

  FILE* cur_file; 

  char file_name[256];

  char i_buffer[4];
  
  for ( i = 0; i< s.n_clusters; i++) {

    file_name[0] = 0;
    strcat(file_name, prefix);
    strcat(file_name, "-");
    sprintf(i_buffer,"%03d", i);
    strcat(file_name, i_buffer);

    cur_file = fopen(file_name, "w");

    for (j = 0; j< s.clusters[i].n_members; j++) {
      fprintf(cur_file,">sequence_%i\n",s.clusters[i].members[j]);
      for ( k = 0;
	    k < (strlen(ds.sequences[s.clusters[i].members[j]])-1);
	    k++) {
	if( k != 0 && k%50 == 0 ) {
	  fprintf(cur_file, "\n");
	  fputc(ds.sequences[s.clusters[i].members[j]][k],cur_file);
	} else {
	  fputc(ds.sequences[s.clusters[i].members[j]][k],cur_file);
	}
      }
      if ( (strlen(ds.sequences[s.clusters[i].members[j]])-1) != 0 &&
	   (strlen(ds.sequences[s.clusters[i].members[j]])-1) %50 == 0 ) {
	fprintf(cur_file, "\n");
	fputc(ds.sequences[s.clusters[i].members[j]][k],cur_file);
      } else {
	fputc(ds.sequences[s.clusters[i].members[j]][k],cur_file);
      }
      fprintf(cur_file, "\n");
    }
    fclose(cur_file);
  }
}
