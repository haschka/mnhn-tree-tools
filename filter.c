#include<stdio.h>
#include<pthread.h>
#include<stdlib.h>
#include<string.h>
#include"smith-waterman.h"
#include"dataset.h"

typedef struct {
  unsigned int max_dist;
  dataset* ds;
  char* target;
  size_t target_length;
  size_t start;
  size_t stride;
  char* is_in_dataset;
} thread_handle_filter_SW;

dataset filter_ds_by_size(dataset in_set, size_t size,
			  size_t plus, size_t minus) {

  int i,counter;

  dataset out_set;

  out_set.sequences = (char**)malloc(sizeof(char*)*in_set.n_values);
  out_set.sequence_lengths = (size_t*)malloc(sizeof(size_t)*in_set.n_values);

  counter = 0;
  for(i=0;i<in_set.n_values;i++) {
    if(in_set.sequence_lengths[i] <= size+plus &&
       in_set.sequence_lengths[i] >= size-minus) {

      out_set.sequences[counter] =
	(char*)malloc(sizeof(char)*in_set.sequence_lengths[i]+1);
      memcpy(out_set.sequences[counter],in_set.sequences[i],
	     sizeof(char)*in_set.sequence_lengths[i]);
      out_set.sequences[counter][in_set.sequence_lengths[i]] = 0;

      out_set.sequence_lengths[counter] = in_set.sequence_lengths[i];
      counter++;
    }
  }
  out_set.n_values = counter;
  out_set.sequences =
    (char**)realloc(out_set.sequences,sizeof(char*)*out_set.n_values);
  out_set.sequence_lengths =
    (size_t*)realloc(out_set.sequence_lengths,sizeof(size_t)*out_set.n_values);

  return(out_set);
}

void* filter_by_SW_distance_thread(void* arg) {

  thread_handle_filter_SW* th = (thread_handle_filter_SW*)arg;

  int i;
  int current_score;
  
  for(i = th->start; i< th->start+th->stride; i++) {

    current_score = score(th->ds->sequences[i],
			  th->target,
			  th->ds->sequence_lengths[i],
			  th->target_length,NULL);
    
    if(th->max_dist >= current_score ) {
      th->is_in_dataset[i] = 1;
    }
  }
}
			     
dataset filter_by_SW_distance(dataset in_set,
			      char* target,
			      size_t target_length,
			      unsigned int max_distance,
			      int n_threads) {

  int i,counter;

  dataset out_set;

  thread_handle_filter_SW* th =
    (thread_handle_filter_SW*)malloc(sizeof(thread_handle_filter_SW)
				     *n_threads);

  pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);

  size_t stride = in_set.n_values/n_threads;
  size_t stride_rest = in_set.n_values%n_threads;

  char* is_in_dataset = (char*)malloc(sizeof(char)*in_set.n_values);

  memset(is_in_dataset,0,sizeof(char)*in_set.n_values);
  
  for(i=0;i<n_threads;i++) {

    th[i].ds = &in_set;
    th[i].target = target;
    th[i].target_length = target_length;
    th[i].start = i*stride;
    th[i].stride = stride;
    if(i == n_threads-1) th[i].stride+=stride_rest;
    th[i].is_in_dataset = is_in_dataset;
    th[i].max_dist = max_distance;
    pthread_create(threads+i, NULL, filter_by_SW_distance_thread, th+i);
  }

  for(i=0;i<n_threads;i++) {
    pthread_join(threads[i],NULL);
  }

  out_set.sequences = (char**)malloc(sizeof(char*)*in_set.n_values);
  out_set.sequence_lengths = (size_t*)malloc(sizeof(size_t)*in_set.n_values);
    
  counter = 0;
  for(i=0;i<in_set.n_values;i++) {
    if(is_in_dataset[i] == 1) {

      out_set.sequences[counter] =
	(char*)malloc(sizeof(char)*in_set.sequence_lengths[i]+1);
      memcpy(out_set.sequences[counter],in_set.sequences[i],
	     sizeof(char)*in_set.sequence_lengths[i]);
      out_set.sequences[counter][in_set.sequence_lengths[i]] = 0;

      out_set.sequence_lengths[counter] = in_set.sequence_lengths[i];
      counter++;
    }
  }
  out_set.n_values = counter;
  out_set.sequences =
    (char**)realloc(out_set.sequences,sizeof(char*)*out_set.n_values);
  out_set.sequence_lengths =
    (size_t*)realloc(out_set.sequence_lengths,sizeof(size_t)*out_set.n_values);

  free(is_in_dataset);
  free(threads);
  free(th);
  
  return(out_set);
}
