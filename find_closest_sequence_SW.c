#include<unistd.h>
#include<stdlib.h>
#include<limits.h>
#include<stdio.h>
#include<pthread.h>
#include"smith-waterman.h"
#include"dataset.h"
#include"cluster.h"


void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

typedef struct{
  int n_threads;
  int thread_num;
  int* partial_scores;
  int* partial_minimal_sequence_indicies;
  dataset* ds;
  dataset* query;
} search_thread_handler;

void* search_thread(void* arg) {

  search_thread_handler* th = (search_thread_handler*)arg;

  int stride; 

  int minimum = INT_MAX;
  int current_score;
  int minimal_sequence_index = -1;
 
  int i,start,end;
  
  if (th->thread_num == 0) {
    stride = th->ds->n_values/th->n_threads+th->ds->n_values%th->n_threads;
    start = 0;
    end = stride;
  } else {
    stride = th->ds->n_values/th->n_threads;
    start = th->ds->n_values%th->n_threads + th->thread_num*stride;
    end = start+stride;
  }

  for(i=start;i<end;i++) {
    current_score = score(th->ds->sequences[i],
			  th->query->sequences[0],
			  th->ds->sequence_lengths[i],
			  th->query->sequence_lengths[0],
			  NULL);

    if (current_score < minimum) {
      minimum = current_score;
      minimal_sequence_index = i;
    }
  }
  th->partial_scores[th->thread_num] = minimum; 
  th->partial_minimal_sequence_indicies[th->thread_num] =
    minimal_sequence_index;
}

int main(int argc, char** argv) {

  int i;
  
  FILE * fasta_f;
  FILE * fasta_q;
  
  int n_threads;
  
  dataset ds;
  dataset query;

  cluster cl;

  char filename[256];

  pthread_t* threads;
  search_thread_handler* th;
  int* partial_score;
  int* partial_indicies;

  int minimum = INT_MAX;
  
  int current_score, current_index, minimal_sequence_index;

  if(argc < 5) {
    printf("Arguments are: \n"
	   "[file] Fasta file of the database to find the closest sequence in\n"
	   "[file] Fasta file containing the query sequence \n"
	   " [int] number of threads to use for this search \n"
	   "[file] output fasta containing the closest sequence in the"
	   " dataset\n");
    return(1);
  }
  
  sscanf(argv[3],"%i",&n_threads);
  sscanf(argv[4],"%s",filename);

  threads = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);
  th = (search_thread_handler*)malloc(sizeof(search_thread_handler)*n_threads);
  partial_score = (int*)malloc(sizeof(int)*n_threads);
  partial_indicies = (int*)malloc(sizeof(int)*n_threads);
  
  cl.members = (int*)malloc(sizeof(int));
  cl.n_members = 1;
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (fasta_q = fopen(argv[2], "r"))) file_error(argv[2]);
  query = dataset_from_fasta(fasta_q);
  fclose(fasta_q);

  for(i=0;i<n_threads;i++) {
    th[i].n_threads = n_threads;
    th[i].thread_num = i;
    th[i].partial_scores = partial_score;
    th[i].partial_minimal_sequence_indicies = partial_indicies;
    th[i].ds = &ds;
    th[i].query = &query;

    pthread_create(threads+i, NULL, search_thread, th+i);
  }

  for(i=0;i<n_threads;i++){
    pthread_join(threads[i],NULL);
  }
 
  for(i=0;i<n_threads;i++) {
    current_score = partial_score[i];
    current_index = partial_indicies[i];

    if (current_score < minimum) {
      minimum = current_score;
      minimal_sequence_index = current_index;
    }
  }

  cl.members[0] = minimal_sequence_index;
  printf("Smith Waterman distance = %i \n", minimum);
  create_single_cluster_file(filename, cl, ds);
}
  
