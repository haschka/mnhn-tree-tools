#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include"dbscan.h"

dataset dataset_from_fasta(FILE* in) {

  char* line = NULL;
  size_t line_size = 0;

  int sequences = 0;

  dataset ds;

  char linebuffer[2000];
  size_t linebuffer_length;
  size_t sequence_length;
  int i;
  
  rewind(in);

  while ( -1 != getline(&line, &line_size, in) ) {
    if( line[0] == '>' ) sequences++;
  }

  rewind(in);

  ds.sequences = (char**)malloc(sizeof(char*)*sequences);
  ds.n_values = sequences;

  sequences = -1;
  
  while ( -1 != getline(&line, &line_size, in) ) {
    if ( line[0] == '>') {
      sequences++;
      ds.sequences[sequences] = (char*)malloc(sizeof(char)*1001);
      ds.sequences[sequences][0] = 0;
      sequence_length = 0;
    } else {
      sscanf(line,"%s",linebuffer);

      linebuffer_length = strlen(linebuffer);
      sequence_length += linebuffer_length;
      if(sequence_length > 1000) {
	ds.sequences[sequences] =
	  (char*)realloc(ds.sequences[sequences],
			 sizeof(char)*(sequence_length+1));
      }
      for(i=0; i < linebuffer_length; i++) {
	linebuffer[i] = toupper(linebuffer[i]);
      }
      strcat(ds.sequences[sequences],linebuffer);
    }
  }
  return(ds);
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
      
int main(int argc, char** argv) {
  
  dataset ds;

  size_t dimensions;
  float epsilon;
  int minpts;
  
  FILE* input = fopen(argv[1], "r");

  split_set set_of_clusters;

  if(argc < 4) {
    printf("Arguments are: \n"
	   " Sequences in FASTA \n"
	   " Epsilon \n"
	   " minPoints \n"
	   " output-files prefix \n" );
    return(1);
  }
  
  sscanf(argv[2], "%f", &epsilon);
  sscanf(argv[3], "%i", &minpts);
  
  ds = dataset_from_fasta(input);

  set_of_clusters = dbscan(ds, epsilon, minpts);

  printf("%u clusters obtained \n", set_of_clusters.n_clusters);

  if (set_of_clusters.n_clusters < 500) {
    create_cluster_files(argv[4], set_of_clusters, ds);
  }
  fclose(input);
}
  
