#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pthread.h> 
#include"dataset.h"

typedef struct {
  size_t id;
  size_t** thread_local_sites;
  size_t* n_local_sites;
  char* sequence;
  dataset * ds;
  size_t start;
  size_t stride;
} thread_handler_digest;


void* digest_thread (void* arg) {

  thread_handler_digest* th = (thread_handler_digest*)arg;
 
  size_t i;

  size_t* n_sites =  th->n_local_sites+(th->id);

  int id = th->id;

  char* sequence = th->sequence;
  
  n_sites[0] = 0;
  
  th->thread_local_sites[id] = (size_t*)malloc(sizeof(size_t)*1000);

  for(i=th->start;i<(th->start+th->stride);i++) {
    if (
#if defined(_digest_XbaI)
	sequence[i-1] == 'T' &&
	sequence[  i] == 'C' &&
	sequence[i+1] == 'T' &&
	sequence[i+2] == 'A' &&
	sequence[i+3] == 'G' &&
	sequence[i+4] == 'A'
#elif defined(_digest_XmnI)
	sequence[i-5] == 'G' &&
	sequence[i-4] == 'A' &&
	sequence[i-3] == 'A' &&
	/*sequence[i-2] == 'N' &&*/
	/*sequence[i-1] == 'N' &&*/
	/*sequence[  i] == 'N' &&*/
	/*sequence[i+1] == 'N' &&*/
	sequence[i+2] == 'T' &&
	sequence[i+3] == 'T' &&
	sequence[i+4] == 'C'
#endif
	) {
      n_sites[0]++;
      if(n_sites[0] %1000 == 0) {
	th->thread_local_sites[id] =
	  (size_t*)realloc(th->thread_local_sites[id],
			   sizeof(size_t)*(n_sites[0]+1000));
      }
      th->thread_local_sites[id][n_sites[0]-1] = i;
    }
  }
}


dataset
#if defined(_digest_XbaI)
digest_XbaI
#elif defined(_digest_XmnI)
digest_XmnI
#endif
(char * sequence, size_t sequence_length, int n_threads) {

  thread_handler_digest* th =
    (thread_handler_digest*)malloc(sizeof(thread_handler_digest)*n_threads);
  pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);
  size_t ** thread_local_sites = (size_t**)malloc(sizeof(size_t*)*n_threads);
  size_t * n_local_sites = (size_t*)malloc(sizeof(size_t*)*n_threads);
  
  size_t i;
  size_t* sites;
  
#if defined(_digest_XbaI)
  size_t bases_prior_cut = 1;
  size_t bases_post_cut = 5;
#elif defined(_digest_XmnI)
  size_t bases_prior_cut = 5;
  size_t bases_post_cut = 5;
#endif
  
  size_t n_sites = 0;
  
  size_t max_size;

  size_t pos;
  
  dataset ds;

  size_t scan_length = (sequence_length-(bases_prior_cut+bases_post_cut)); 
  
  size_t stride = scan_length/n_threads;
  size_t stride_rest = scan_length%n_threads;

  for(i=0;i<n_threads;i++) {
    th[i].id = i;
    th[i].n_local_sites = n_local_sites;
    th[i].thread_local_sites = thread_local_sites;
    th[i].sequence= sequence;
    th[i].start=bases_prior_cut+i*stride;
    th[i].stride=stride;
    if(i == n_threads - 1) th[i].stride += stride_rest;

    pthread_create(threads+i,NULL,digest_thread,th+i);
  }

  for(i=0;i<n_threads;i++) {
    pthread_join(threads[i],NULL);
  }

  free(th);
  free(threads);
  
  for(i=0;i<n_threads;i++) {
    n_sites+=n_local_sites[i];
  }
  
  sites = (size_t*)malloc(sizeof(size_t)*n_sites);
  pos = 0;
  for(i=0;i<n_threads;i++) {
    if(n_local_sites[i] > 0) {
      memcpy(sites+pos,thread_local_sites[i],sizeof(size_t)*n_local_sites[i]);
      pos += n_local_sites[i];
    }
    free(thread_local_sites[i]);
  }
  free(n_local_sites);
  free(thread_local_sites);
    
  ds.n_values = n_sites+1;
  ds.sequences = (char**)malloc(sizeof(char*)*(ds.n_values));
  ds.sequence_lengths = (size_t*)malloc(sizeof(size_t)*(ds.n_values));
    
  
  if(n_sites == 0) {
    ds.sequences[0] = (char*)malloc(sizeof(char)*(sequence_length+1));
    memcpy(ds.sequences[0],sequence,sequence_length);
    ds.sequences[0][sequence_length] = 0;
    ds.sequence_lengths[0] = sequence_length;
  }

  if(n_sites > 0) {
    
    /*first*/
    ds.sequences[0] = (char*)malloc(sizeof(char)*(sites[0])+1);
    memcpy(ds.sequences[0],sequence,sizeof(char)*(sites[0]));
    ds.sequences[0][sites[0]] = 0;
    ds.sequence_lengths[0] = sites[0];
    /*last*/
    ds.sequence_lengths[n_sites] = sequence_length-sites[n_sites-1];
    ds.sequences[n_sites] =
      (char*)malloc(sizeof(char)*(ds.sequence_lengths[n_sites]+1));
    memcpy(ds.sequences[n_sites],sequence+sites[n_sites-1],sizeof(char)*
	   ds.sequence_lengths[n_sites]);
    ds.sequences[n_sites][ds.sequence_lengths[n_sites]] = 0; 
    /*inbetween*/
    for(i=0;i<n_sites-1;i++) {
      ds.sequence_lengths[i+1] = sites[i+1]-sites[i];
      ds.sequences[i+1] = (char*)malloc(sizeof(char)*
					ds.sequence_lengths[i+1]+1);
      
      memcpy(ds.sequences[i+1],sequence+sites[i],ds.sequence_lengths[i+1]);
      ds.sequences[i+1][ds.sequence_lengths[i+1]] = 0;
    }
  }

  max_size = 0;
  for(i=0;i<ds.n_values;i++) {
    if(ds.sequence_lengths[i] > max_size) max_size = ds.sequence_lengths[i];
  }
  ds.max_sequence_length = max_size;

  free(sites);
  
  return(ds);
}
	
    
