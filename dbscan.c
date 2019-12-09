#include<stdlib.h>
#include<string.h>
#include<stdio.h>

#if defined(_SCAN_L1)
#include<math.h>
#endif

#include"binary_array.h"
#include"dataset.h"
#include"cluster.h"
#include"dbscan.h"

#if defined(_SCAN_SMITH_WATERMAN)
#include"smith-waterman.h"
#endif

#if defined(__AVX__)
#include<immintrin.h>
#elif defined(__SSE__)
#include<xmmintrin.h>
#endif

static inline int compare(const void* a , const void* b) {
  int na = *(int*)a;
  int nb = *(int*)b;
  return (na-nb);
}

#if defined(_SCAN_SMITH_WATERMAN)
static inline neighbors region_query(int point, float epsilon, dataset ds) {

  int i,j,k;

  neighbors nb;

  float* a_point;
  float* b_point;

  float distance;
  float a_minus_b;

  int* work;

  nb.members = (int*)malloc(sizeof(int)*ds.n_values);
  nb.n_members = 0;

  work =
    (int*)malloc(sizeof(int)
		 *(ds.max_sequence_length+1)*(ds.max_sequence_length+1));
  
  for(i=0;i<ds.n_values;i++) {

    distance = (float)score(ds.sequences[i],ds.sequences[point],
			    ds.sequence_lengths[i],ds.sequence_lengths[point],
			    work);

    if(distance <= epsilon) {
      nb.members[nb.n_members]=i;
      nb.n_members++;
    }
  }

  /*  nb.members = (int*)realloc((void*)nb.members,sizeof(int)*nb.n_members); */
  qsort(nb.members,nb.n_members,sizeof(int),compare);
  free(work);
  return(nb);
}

#elif defined(_SCAN_L2) || defined(_SCAN_L1)

#if defined(_SCAN_L1)
static int abs_sign_i = 0x7fffffff;
#endif 
static inline neighbors region_query(int point, float epsilon, dataset ds) {

  int i,j,k;

  neighbors nb;

  float* a_point;
  float* b_point;

  float epsilon_square = epsilon*epsilon;

  float distance;
  float a_minus_b;

  a_point = (float*)malloc(sizeof(float)*ds.n_dimensions);
  b_point = (float*)malloc(sizeof(float)*ds.n_dimensions);

#if defined(__AVX__)

  __m256* point_v;
  __m256 current_v;
  __m256 a_minus_b_v;
  __m256 dist_squared;
  __m256 cmp_v;

  int i_eight;
  
  float result[8] __attribute__((aligned(32)));

  int end;
#if defined(_SCAN_L1)
  __m256 abs_sign;
  float* abs_sign_f = (float*)&abs_sign_i;
  abs_sign = _mm256_set1_ps(abs_sign_f[0]);
#endif
  posix_memalign((void**)&point_v,32,32*ds.n_dimensions);
  
  
#elif defined(__SSE__)

  __m128* point_v;
  __m128 current_v;
  __m128 a_minus_b_v;
  __m128 dist_squared;
  __m128 cmp_v;

  int i_four;

  float result[4] __attribute__((aligned(16)));

  int end;
#if defined(_SCAN_L1)
  __m128 abs_sign;
  float* abs_sign_f = (float*)&abs_sign_i;
  abs_sign = _mm_set1_ps(abs_sign_f[0]);
#endif
  posix_memalign((void**)&point_v,16,16*ds.n_dimensions);
  
#endif

  nb.members = (int*)malloc(sizeof(int)*ds.n_values);
  nb.n_members = 0;

#if defined(__AVX__)

  for(k = 0;k < ds.n_dimensions; k++) {
    point_v[k] = _mm256_set1_ps(ds.values[k][point]);
  }
  
  for(i=0;i<ds.n_values/8+(ds.n_values%8 != 0);i++) {

    i_eight = i*8;

    dist_squared = _mm256_setzero_ps();

    for(k = 0; k < ds.n_dimensions; k++) {
      
      current_v = _mm256_load_ps(ds.values[k]+i_eight);
      a_minus_b_v = _mm256_sub_ps(point_v[k], current_v);

#if defined(_SCAN_L2)
#ifdef __FMA__

      dist_squared = _mm256_fmadd_ps(a_minus_b_v,
				     a_minus_b_v,
				     dist_squared);

#else
      dist_squared = _mm256_add_ps(_mm256_mul_ps(a_minus_b_v,
						 a_minus_b_v),
				   dist_squared);
#endif
#elif defined(_SCAN_L1)
       dist_squared = _mm256_add_ps(_mm256_and_ps(a_minus_b_v,
						  abs_sign),
				    dist_squared);
#endif
    }
    _mm256_store_ps(result,dist_squared);

    if (__builtin_expect(i != ds.n_values/8,1)) {
      end = 8;
    } else {
      end = ds.n_values%8;
    }

    for(j=0;j<end;j++) {
      if (
#if defined(_SCAN_L2)
	  result[j] <= epsilon_square 
#elif defined(_SCAN_L1)
	  result[j] <= epsilon 
#endif
	  ) {
	nb.members[nb.n_members]=i_eight+j;
	nb.n_members++;
      }
    }
  }

#elif defined(__SSE__)

  for(k = 0; k < ds.n_dimensions; k++) {
    point_v[k] = _mm_set1_ps(ds.values[k][point]);
  }

  for(i=0;i<ds.n_values/4+(ds.n_values%4 != 0);i++) {

    i_four = i*4;

    dist_squared = _mm_setzero_ps();

    for(k = 0; k < ds.n_dimensions; k ++) {
      current_v = _mm_load_ps(ds.values[k]+i_four);
      
      a_minus_b_v = _mm_sub_ps(point_v[k],current_v);

#if defined(_SCAN_L2) 
      dist_squared = _mm_add_ps(_mm_mul_ps(a_minus_b_v,a_minus_b_v),
				dist_squared);
#elif defined(_SCAN_L1)
      dist_squared = _mm_add_ps(_mm_and_ps(a_minus_b_v,
				   abs_sign),
			dist_squared);
#endif
    }

    _mm_store_ps(result,dist_squared);

    if (__builtin_expect(i != ds.n_values/4,1)) {
      end = 4;
    } else {
      end = ds.n_values%4;
    }

    for(j=0;j<end;j++) {
      if(
#if defined (_SCAN_L2)
	 result[j] <= epsilon_square  
#elif defined (_SCAN_L1)
	 result[j] <= epsilon
#endif
	 ) {
	nb.members[nb.n_members]=i_four+j;
	nb.n_members++;
      }
    }
  }
#else
  for(k = 0; k< ds.n_dimensions; k++) {
    a_point[k] = ds.values[k][point];
  }

  nb.n_members = 0;
  for(i=0;i<ds.n_values;i++) {

    for(k = 0; k < ds.n_dimensions; k++) {
      b_point[k] = ds.values[k][i];
    }
    distance = 0;

    for(j=0;j< ds.n_dimensions;j++) {
      a_minus_b=(a_point[j]-b_point[j]);
#if defined (_SCAN_L2)
      distance += a_minus_b*a_minus_b;
#elif defined (_SCAN_L1)
      distance += fabsf(a_minus_b);
#endif
    }
    if(
#if defined (_SCAN_L2)
       distance <= epsilon_square
#elif defined(_SCAN_L1)
       distance <= epsilon
#endif
       ) {
      nb.members[nb.n_members]=i;
      nb.n_members++;
    }
  }
#endif

#if defined(__AVX__) || defined (__SSE__)
  free(point_v);
#endif
  free(a_point);
  free(b_point);
  
  /*  nb.members = (int*)realloc((void*)nb.members,sizeof(int)*nb.n_members); */
  qsort(nb.members,nb.n_members,sizeof(int),compare);
  return(nb);
}
#endif


static inline void expand_cluster(int point,
				  neighbors nb,
				  cluster* cl,
				  float epsilon,
				  int minpts,
				  char* visited,
				  int* cluster_member,
				  dataset ds) {

  int i,j;
  neighbors nb_of_nb;
  int * merge;
  int * right_in_left;
  int n_right_in_left;
  int merge_counter;
  int current_point;

  int add_to_unsorted_counter;
  int right_in_left_counter;

  int max_members_nb_nb_of_nb;
  int left;
  int right;
  int n_left;
  int n_right;

  neighbors nb_unsorted;

  nb_unsorted.members = (int*)malloc(sizeof(int)*ds.n_values);

  merge = (int*)malloc(ds.n_values*sizeof(int));

  memcpy(nb_unsorted.members,nb.members,sizeof(int)*nb.n_members);
  nb_unsorted.n_members = nb.n_members;

  cl->members = (int*)malloc(sizeof(int)*ds.n_values);
  cl->members[cl->n_members] = point;
  cluster_member[point] = 1;
  cl->n_members++;

  for(i=0;i < nb_unsorted.n_members;i++) {
    if(!get_value_in_binary_array_at_index(visited,nb_unsorted.members[i])) {
      set_value_in_binary_array_at_index(visited,nb_unsorted.members[i]);
      nb_of_nb = region_query(nb_unsorted.members[i], epsilon, ds);
      if (nb_of_nb.n_members >= minpts) {
	/*merge = (int*)malloc((nb.n_members+nb_of_nb.n_members)
	 *sizeof(int));*/

	right_in_left = (int*)malloc(nb_of_nb.n_members*sizeof(int));
	n_left = n_right = merge_counter = 0;
	n_right_in_left = 0;
	while(n_left < nb.n_members && n_right < nb_of_nb.n_members) {
	  left = nb.members[n_left];
	  right = nb_of_nb.members[n_right];
	  if( left < right ) {
	    merge[merge_counter] = left;
	    merge_counter++;
	    n_left++;
	  }
	  if( right < left ) {
	    merge[merge_counter] = right;
	    merge_counter++;
	    n_right++;
	  }
	  if( right == left ) {
	    merge[merge_counter] = left;
	    n_left++;
	    n_right++;
	    merge_counter++;
	    right_in_left[n_right_in_left] = left;
	    n_right_in_left++;
	  }
	}
	if(n_left != nb.n_members || n_right != nb_of_nb.n_members) {
	  if(n_left == nb.n_members) { /* Hit the left wall */
	    for(j=n_right;j<nb_of_nb.n_members;j++) {
	      merge[merge_counter] = nb_of_nb.members[j];
	      merge_counter++;
	    }
	  }
	  if(n_right == nb_of_nb.n_members) { /*Hit the right wall */
	    for(j=n_left;j<nb.n_members;j++) {
	      merge[merge_counter] = nb.members[j];
	      merge_counter++;
	    }
	  }
	}

	/*	nb_unsorted.members = (int*)realloc((void*)nb_unsorted.members,
		merge_counter*sizeof(int));*/

	add_to_unsorted_counter = 0;
	right_in_left_counter = 0;
	/* if there are rights in left */
	if (n_right_in_left) {
	  for(j=0;j<nb_of_nb.n_members;j++) {

	    if(right_in_left[right_in_left_counter] != nb_of_nb.members[j]) {

	      nb_unsorted.members[nb_unsorted.n_members
				  +add_to_unsorted_counter] =
		nb_of_nb.members[j];

	      add_to_unsorted_counter++;

	    } else {
	      right_in_left_counter++;
	    }
	  }
	  /* if there are no rights in left */
	} else {
	  for(j=0;j<nb_of_nb.n_members;j++) {
	    nb_unsorted.members[nb_unsorted.n_members+j] = nb_of_nb.members[j];
	  }
	}
	free(right_in_left);
	free(nb_of_nb.members);
	nb_unsorted.n_members = merge_counter;
	/*merge = (int*)realloc((void*)merge,merge_counter*sizeof(int));*/
	memcpy(nb.members,merge,sizeof(int)*merge_counter);
	//nb.members = merge;

	nb.n_members = merge_counter;
	/*
	  nb.members = (int*)realloc((void*)nb.members,
	  merge_counter*sizeof(int));
	*/
      }
    }
    if(!cluster_member[nb_unsorted.members[i]]) {
      /*      cl->members = (int*)realloc((void*)(cl->members)
	      ,sizeof(int)*(cl->n_members+1));*/
      cl->members[cl->n_members] = nb_unsorted.members[i];
      cl->n_members++;
      cluster_member[nb_unsorted.members[i]] = 1;

    }
  }
  cl->members=(int*)realloc((void*)(cl->members),sizeof(int)*(cl->n_members));

  free(nb_unsorted.members);
  free(nb.members);
}

split_set
#if defined (_SCAN_L1)
dbscan_L1
#elif defined (_SCAN_L2)
dbscan_L2
#elif defined (_SCAN_SMITH_WATERMAN)
dbscan_SW
#endif
(dataset ds, float epsilon, int minpts) {

  int i;
  neighbors nb;
  split_set ret_val;
  cluster* clusters = (cluster*)malloc(sizeof(cluster)*ds.n_values);
  char* visited = alloc_and_set_zero_binary_array(ds.n_values);
  int n_clusters = 0;
  int* cluster_member = (int*)malloc(sizeof(int)*ds.n_values);
  
  if (epsilon == -1) {
    printf("Warning @dbscan: epsilon value is: %f\n", epsilon);
    epsilon = 0.015;
    printf("Warning @dbscan: epsilon was changed to: %f\n", epsilon);
  }
  if (minpts == -1) {
    printf("Warning @dbscan: minpts value is: %i\n", minpts);
    minpts = 5;
    printf("Warning @dbscan: minpts was changed to: %i\n", minpts);
  }

  ret_val.clusters = (cluster*)malloc(sizeof(cluster)*ds.n_values);

  memset(cluster_member,0,sizeof(int)*ds.n_values);

  for(i=0;i<ds.n_values;i++) {
    //    clusters[i].members = NULL;
    clusters[i].n_members = 0;
    clusters[i].id = i;
  }

  for(i=0;i<ds.n_values;i++) {
    if(!get_value_in_binary_array_at_index(visited,i)) {
      set_value_in_binary_array_at_index(visited,i);
      nb = region_query(i, epsilon, ds);
      if ( minpts < nb.n_members ) {
	expand_cluster(i,nb,clusters+(n_clusters),epsilon,minpts,
		       visited,cluster_member,ds);
	n_clusters++;
      } else {
	free(nb.members);
      }
    }
  }
  //clusters = (cluster*)realloc(clusters,sizeof(cluster)*n_clusters);

  ret_val.n_clusters = 0;
  for(i=0;i<n_clusters;i++) {
    //if(clusters[i].n_members) {
      memcpy(ret_val.clusters+ret_val.n_clusters,clusters+i,sizeof(cluster));
      ret_val.n_clusters++;
      //}
  }
  free(clusters);
  free(visited);
  free(cluster_member);

  return(ret_val);
}


