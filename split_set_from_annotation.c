#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#include"dataset.h"
#include"cluster.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

void print_arguments() {

  printf("Arguments are: \n");
  printf(" [FILE-in]  Fasta file to cluster according to annotations \n");
  printf(" [FILE-in]  annotation file - one string per annotation \n");
  printf("            lines only containing '-' are skipped - outliers \n");
  printf(" [FILE-in]  unique annotations - one string per annotation in \n");
  printf("            split-set order \n");
  printf(" [FILE-out] split-set file to be written \n");
  _exit(1);
}

int main(int argc, char** argv) {

  int i, j;
  
  dataset ds;

  FILE* fasta_f;
  FILE* annotation_f;
  FILE* uniq_annotations_f;
  
  int n_annotations = argc-4;

  char** uniq_annotations;
  char** annotation_table;

  char scan_buffer[2000];
  
  char* current_annotation;
  size_t current_length;
  size_t max_length;
  
  split_set s;

  if(argc < 4) print_arguments();
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (annotation_f = fopen(argv[2], "r"))) file_error(argv[2]);
  if ( NULL == (uniq_annotations_f = fopen(argv[3], "r"))) file_error(argv[3]);

  n_annotations = 0;
  annotation_table = NULL;

  while( 1 == fscanf(uniq_annotations_f,"%s", scan_buffer)) {
    current_annotation = scan_buffer;
    current_length = strlen(current_annotation);
    if (current_length == 1 && current_annotation[0] == '-') {
      continue;
    }
    current_length++;

    annotation_table = (char**)realloc(annotation_table,
				       sizeof(char*)*(n_annotations+1));
    annotation_table[n_annotations] =
      (char*)malloc(sizeof(char)*(current_length));

    memcpy(annotation_table[n_annotations],
	   current_annotation,
	   sizeof(char)*current_length);
    n_annotations++;

  }
  
  s.n_clusters = n_annotations;
  s.clusters = (cluster*)malloc(sizeof(cluster)*n_annotations);

  max_length = 0;
  for(i=0;i<n_annotations;i++) {

    s.clusters[i].n_members = 0;
    s.clusters[i].members = malloc(sizeof(int)*1000);

  }   

  for(i=0;i<ds.n_values;i++) {
    if (1 != fscanf(annotation_f,"%s",scan_buffer)) {
      printf("Annotations file corrupt ! \n");
      return(1);
    }
    if( strlen(scan_buffer) == 1 && scan_buffer[0] == '-' ) {
      continue;
    }
    for(j=0;j<n_annotations;j++) {
      if (!strcmp(scan_buffer,annotation_table[j])) {
	s.clusters[j].n_members++;
	if(s.clusters[j].n_members%1000 == 0) {
	  s.clusters[j].members = realloc(s.clusters[j].members,
					  sizeof(int)
					  *(s.clusters[j].n_members+1000));
	}
	s.clusters[j].members[s.clusters[j].n_members-1] = i;
      }
    }
  }

  for(i=0;i<n_annotations;i++) free(annotation_table[i]);
  free(annotation_table);
  store_split_set(argv[4],s);
  free_split_set_and_associated_clusters(s);
  free_sequences_from_dataset(ds);
  return(0);
}
	
    
    
  
