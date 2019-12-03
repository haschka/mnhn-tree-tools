#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<limits.h>
#include<pthread.h>
#include"dataset.h"
#include"smith-waterman.h"

typedef struct {
  dataset* ds_one;
  dataset* ds_two;
  size_t id;
  size_t* incrementor;
  double* local_sums;
  double* mean;
} thread_handle;

pthread_mutex_t lock;

void* silhouette_from_smith_waterman_datasets_thread(void* th) {

  double internal;
  double external;

  thread_handle* handle = (thread_handle*)th;
  
  dataset* ds_one = handle->ds_one;
  dataset* ds_two = handle->ds_two;
  
  double a,b,s, s_mean;

  size_t i,j;

  int* work;
  int worksize;

  if(ds_one->max_sequence_length >= ds_two->max_sequence_length) {
    worksize = ds_one->max_sequence_length*ds_one->max_sequence_length;
  } else {
    worksize = ds_two->max_sequence_length*ds_two->max_sequence_length;
  }
  work = (int*)malloc(sizeof(int)*worksize);
  
  s_mean = 0;

 thread_continuition:

  pthread_mutex_lock(&lock);
  i = handle->incrementor[0];
  handle->incrementor[0]++;
  pthread_mutex_unlock(&lock);

  if(i < ds_one->n_values) {
    internal = 0.;
    external = 0.;
    for(j=0;j<ds_one->n_values;j++) {
      internal += (double)score(ds_one->sequences[i],ds_one->sequences[j],
				ds_one->sequence_lengths[i],
				ds_one->sequence_lengths[j],work);
    }
    for(j=0;j<ds_two->n_values;j++) {
      external += (double)score(ds_one->sequences[i],ds_two->sequences[j],
				ds_one->sequence_lengths[i],
				ds_two->sequence_lengths[j],work);
    }
    
    a=(internal/((double)ds_one->n_values-1.));
    b=(external/((double)ds_two->n_values));

    if(a<b) {
      s = 1-a/b;
    } else if(a==b) {
      s = 0;
    } else if(a>b) {
      s = b/a-1;
    }
    s_mean += s;
    goto thread_continuition;
  }
  free(work);
  handle->local_sums[handle->id] = s_mean;
}
 
unsigned long* smith_waterman_distances_matrix(dataset ds_one,
						dataset ds_two) {

  unsigned long* matrix =
    (unsigned long*)malloc(sizeof(unsigned long)
			   *ds_one.n_values*ds_two.n_values);

  size_t i,j;
  size_t worksize;
  int* work;
  
  if(ds_one.max_sequence_length >= ds_two.max_sequence_length) {
    worksize = ds_one.max_sequence_length*ds_one.max_sequence_length;
  } else {
    worksize = ds_two.max_sequence_length*ds_two.max_sequence_length;
  }
  work = (int*)malloc(sizeof(int)*worksize);
  
  for(i=0;i<ds_one.n_values;i++) {
    for(j=0;j<ds_two.n_values;j++) {

      matrix[i*ds_two.n_values+j] =
	score(ds_one.sequences[i],ds_two.sequences[j],
	      ds_one.sequence_lengths[i],ds_two.sequence_lengths[j], work);
    }
  }
  free(work);
  return(matrix);
}

void* mean_from_smith_waterman_datasets_thread(void* th) {

  thread_handle* handle = (thread_handle*)th;
  
  dataset* ds_one = handle->ds_one;
  dataset* ds_two = handle->ds_two;

  int sum = 0;

  int* work;
  int worksize;

  size_t i,j;

  if(ds_one->max_sequence_length >= ds_two->max_sequence_length) {
    worksize = ds_one->max_sequence_length*ds_one->max_sequence_length;
  } else {
    worksize = ds_two->max_sequence_length*ds_two->max_sequence_length;
  }
  work = (int*)malloc(sizeof(int)*worksize);

 thread_continuition:

  pthread_mutex_lock(&lock);
  i = handle->incrementor[0];
  handle->incrementor[0]++;
  pthread_mutex_unlock(&lock);

  if(i<ds_one->n_values) {
    for(j=0;j<ds_two->n_values;j++) {

      sum+=
	score(ds_one->sequences[i],ds_two->sequences[j],
	      ds_one->sequence_lengths[i], ds_two->sequence_lengths[j], work);
    }
    goto thread_continuition;
  }
  handle->local_sums[handle->id] = sum;
  free(work);
}


double mean_from_smith_waterman_datasets(dataset ds_one, dataset ds_two,
					 size_t n_threads) {

  size_t i;

  double sum = 0.;
  
  size_t incrementor = 0;

  pthread_t* thread = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);

  thread_handle* th = (thread_handle*)malloc(sizeof(thread_handle)*n_threads);
  double* local_sums = (double*)malloc(sizeof(double)*n_threads);

   for(i=0;i<n_threads;i++) {
    th[i].ds_one = &ds_one;
    th[i].ds_two = &ds_two;
    th[i].id = i;
    th[i].incrementor = &incrementor;
    th[i].local_sums = local_sums;
    
    pthread_create(thread+i, NULL,
		    mean_from_smith_waterman_datasets_thread, th+i);
  }

  for(i=0;i<n_threads;i++) {
    pthread_join(thread[i], NULL);
  }
  for(i=0;i<n_threads;i++) {
    sum += local_sums[i];
  }
  
  free(thread);
  free(th);
  free(local_sums);
  
  return(sum/((double)ds_one.n_values*(double)ds_two.n_values));
}

void* sigma_from_smith_waterman_datasets_thread(void* th) {

  thread_handle* handle = (thread_handle*)th;
  
  dataset* ds_one = handle->ds_one;
  dataset* ds_two = handle->ds_two;

  double mean = handle->mean[0];
  
  int sum = 0;
  double v;

  int* work;
  int worksize;

  size_t i,j;

  if(ds_one->max_sequence_length >= ds_two->max_sequence_length) {
    worksize = ds_one->max_sequence_length*ds_one->max_sequence_length;
  } else {
    worksize = ds_two->max_sequence_length*ds_two->max_sequence_length;
  }
  work = (int*)malloc(sizeof(int)*worksize);

  
 thread_continuition:

  pthread_mutex_lock(&lock);
  i = handle->incrementor[0];
  handle->incrementor[0]++;
  pthread_mutex_unlock(&lock);

  if(i<ds_one->n_values) {
    for(j=0;j<ds_two->n_values;j++) {
      v = score(ds_one->sequences[i],ds_two->sequences[j],
		ds_one->sequence_lengths[i],
		ds_two->sequence_lengths[j], work);
      sum+=((double)v-(double)mean)*((double)v-(double)mean);
    }
    goto thread_continuition;
  }
  handle->local_sums[handle->id] = sum;
  free(work);
}


double sigma_from_smith_waterman_datasets(double mean,
					  dataset ds_one,
					  dataset ds_two,
					  size_t n_threads) {

  size_t i;
  double sigma = 0.;

  size_t incrementor = 0;

  pthread_t* thread = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);

  thread_handle* th = (thread_handle*)malloc(sizeof(thread_handle)*n_threads);
  double* local_sums = (double*)malloc(sizeof(double)*n_threads);

   for(i=0;i<n_threads;i++) {
    th[i].ds_one = &ds_one;
    th[i].ds_two = &ds_two;
    th[i].id = i;
    th[i].mean = &mean;
    th[i].incrementor = &incrementor;
    th[i].local_sums = local_sums;
    
    pthread_create(thread+i, NULL,
		   sigma_from_smith_waterman_datasets_thread, th+i);
  }

  for(i=0;i<n_threads;i++) {
    pthread_join(thread[i], NULL);
  }
  for(i=0;i<n_threads;i++) {
    sigma += local_sums[i];
  }
  
  free(thread);
  free(th);
  free(local_sums);

  return(sqrt(sigma)/(ds_one.n_values+ds_two.n_values));
}

double silhouette_from_smith_waterman_datasets(dataset ds_one,
					       dataset ds_two,
					       size_t n_threads) {

  double internal;
  double external;

  double a,b,s, s_mean = 0;
  
  size_t i,j;

  size_t incrementor = 0;

  pthread_t* thread = (pthread_t*)malloc(sizeof(pthread_t)*n_threads);

  thread_handle* th = (thread_handle*)malloc(sizeof(thread_handle)*n_threads);
  double* local_sums = (double*)malloc(sizeof(double)*n_threads);
  
  for(i=0;i<n_threads;i++) {
    th[i].ds_one = &ds_one;
    th[i].ds_two = &ds_two;
    th[i].id = i;
    th[i].incrementor = &incrementor;
    th[i].local_sums = local_sums;
    
    pthread_create(thread+i, NULL,
		   silhouette_from_smith_waterman_datasets_thread, th+i);
  }

  for(i=0;i<n_threads;i++) {
    pthread_join(thread[i], NULL);
  }
  for(i=0;i<n_threads;i++) {
    s_mean += local_sums[i];
  }

  free(thread);
  free(th);
  free(local_sums);
  
  return(s_mean/(double)ds_one.n_values);
}

double mean_from_smith_waterman_distance_matrix(unsigned long* matrix,
						dataset ds_one,
						dataset ds_two) {

  size_t i;
  
  double sum = 0.;
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    sum+=(double)matrix[i];
  }
  return(sum/((double)ds_one.n_values*(double)ds_two.n_values));
  
}

double sigma_from_smith_waterman_distance_matrix(unsigned long* matrix,
						 double mean,
						 dataset ds_one,
						 dataset ds_two) {

  size_t i;
  double sigma = 0.;
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    sigma+=((double)matrix[i]-(double)mean)*((double)matrix[i]-(double)mean);
  }
  return(sqrt(sigma)/((double)ds_one.n_values+(double)ds_two.n_values));
}

void shortest_longest_distance_in_matrix(unsigned long* matrix,
					 unsigned long* min,
					 unsigned long* max,
					 dataset ds_one, dataset ds_two) {

  min[0] = ULONG_MAX;
  max[0] = 0;

  size_t i;
  
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    if(matrix[i] < min[0]) min[0] = matrix[i];
    if(matrix[i] > max[0]) max[0] = matrix[i];
  }
}
				       
