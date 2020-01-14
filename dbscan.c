#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<unistd.h>

#if defined(_SCAN_SMITH_WATERMAN_GPU)
#ifdef __APPLE__
#include<OpenCL/OpenCL.h>
#else
#include<CL/opencl.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#endif

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

#if defined(_SCAN_SMITH_WATERMAN_GPU)

char* load_program_source(const char* filename) {
  	
  struct stat statbuf;
  FILE *fh; 
  char *source; 
  
  fh = fopen(filename, "r");
  if (fh == 0)
    return 0; 
  
  stat(filename, &statbuf);
  source = (char *) malloc(statbuf.st_size + 1);
  fread(source, statbuf.st_size, 1, fh);
  source[statbuf.st_size] = '\0'; 
  
  return source; 
}

opencl_stuff opencl_initialization(dataset ds) {

  int i;
  
  opencl_stuff ocl;

  cl_platform_id* platform; 
  cl_context_properties* contprop;
  cl_program* programs;
  
  cl_uint num_gpus;
  cl_uint num_acc;
  cl_uint num_cpus;

  int num_devs;

  cl_int err;
  
  cl_device_id* devices;

  char BuildErrorLog[2048];
  size_t BuildErrorLength;

  char clSourceFile[] = "smith-waterman-gpu.cl";
  char* clSource;
  
  size_t acc_distance_size;
  
  int seq_offset = 300;

  char* gpu_sequence_dataset_buffer=
    (char*)malloc(sizeof(char)*300*ds.n_values);

  for(i=0;i<ds.n_values;i++) {
    memcpy(gpu_sequence_dataset_buffer+300*i,
	   ds.sequences[i],ds.sequence_lengths[i]+1);
  }

  platform = (cl_platform_id*)malloc(sizeof(cl_platform_id));
  
  clGetPlatformIDs(1, platform,NULL);
  err = clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_gpus);
  if (err == CL_DEVICE_NOT_FOUND) num_gpus = 0;
  err = clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_CPU, 0, NULL, &num_cpus);
  if (err == CL_DEVICE_NOT_FOUND) num_cpus = 0;
  err = clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_ACCELERATOR, 0, NULL, &num_acc);
  if (err == CL_DEVICE_NOT_FOUND) num_acc = 0;

  if (num_acc == 0 && num_gpus == 0 && num_cpus == 0) {
    printf("No opencl devices found");
    _exit(1);
  }
  
  if ( num_gpus > 0 ) {
    devices = (cl_device_id*)malloc(sizeof(cl_device_id)*num_gpus);
    clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_GPU, num_gpus, devices, NULL);
    num_devs = num_gpus;
  } else if (num_acc > 0) {
    devices = (cl_device_id*)malloc(sizeof(cl_device_id)*num_acc);
    clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_ACCELERATOR,
		   num_acc, devices, NULL);
    num_devs = num_acc;
  } else if (num_cpus > 0 ) {
    devices = (cl_device_id*)malloc(sizeof(cl_device_id)*num_cpus);
    clGetDeviceIDs(platform[0], CL_DEVICE_TYPE_CPU, num_cpus, devices, NULL);
    num_devs = num_cpus;
  }

  ocl.cmdq = (cl_command_queue*)malloc(sizeof(cl_command_queue)*num_devs);
  ocl.contexts = (cl_context*)malloc(sizeof(cl_context)*num_devs);
  ocl.kernel = (cl_kernel*)malloc(sizeof(cl_kernel)*num_devs);
  programs = (cl_program*)malloc(sizeof(cl_program)*num_devs);

  ocl.acc_sequences = (cl_mem*)malloc(sizeof(cl_mem)*num_devs);
  ocl.acc_sequence_lengths = (cl_mem*)malloc(sizeof(cl_mem)*num_devs);
  ocl.acc_distances = (cl_mem*)malloc(sizeof(cl_mem)*num_devs);

  ocl.local_distances = (int*)malloc(sizeof(int)*ds.n_values);

  contprop = (cl_context_properties*)malloc(sizeof(cl_context_properties)*3);
  
  contprop[0] = CL_CONTEXT_PLATFORM;
  contprop[1] = (cl_context_properties)platform[0];
  contprop[2] = 0;

  clSource = load_program_source(clSourceFile);
  
  for(i=0;i<num_devs;i++) {
    ocl.contexts[i] = clCreateContext(contprop, 1, devices+i, NULL, NULL, &err);
    if(err != CL_SUCCESS) {
      printf("Context Creation failed! OCL-ERROR-CODE: %i\n", err);
    }
    ocl.cmdq[i] = clCreateCommandQueue(ocl.contexts[i], devices[i], 0, NULL);

    programs[i] = clCreateProgramWithSource(ocl.contexts[i],1,
					    (const char**)&clSource,NULL,
					    &err);
    if( err != CL_SUCCESS) {
      printf("Program creation failed! \n");
    }
    err = clBuildProgram(programs[i],0,NULL, NULL, NULL, NULL);

    if( err != CL_SUCCESS) {
      clGetProgramBuildInfo(programs[i],devices[i],CL_PROGRAM_BUILD_LOG,
			    2048,BuildErrorLog,&BuildErrorLength);

      printf("Build Error Log: \n%s\n",BuildErrorLog); 
    }

    //free(clSource);
    
    ocl.kernel[i] = clCreateKernel(programs[i],"gpuwaterman",&err);

    ocl.acc_sequence_lengths[i] = clCreateBuffer(ocl.contexts[i],
					     CL_MEM_READ_WRITE,
					     sizeof(size_t)*ds.n_values,
					     NULL,&err);
    err = clEnqueueWriteBuffer(ocl.cmdq[i],
			       ocl.acc_sequence_lengths[i],
			       CL_TRUE,
			       0,
			       sizeof(size_t)*ds.n_values,
			       ds.sequence_lengths,
			       0,NULL,NULL);
    
    ocl.acc_sequences[i] = clCreateBuffer(ocl.contexts[i],
				      CL_MEM_READ_WRITE,
				      sizeof(char)*300*ds.n_values,
				      NULL,&err);

    err = clEnqueueWriteBuffer(ocl.cmdq[i],
			       ocl.acc_sequences[i],
			       CL_TRUE,
			       0,
			       sizeof(char)*300*ds.n_values,
			       gpu_sequence_dataset_buffer,
			       0,NULL,NULL);

    if (i == 0) {
      acc_distance_size = ds.n_values/num_devs + ds.n_values%num_devs;
    } else {
      acc_distance_size = ds.n_values/num_devs;
    }
    ocl.acc_distances[i] = clCreateBuffer(ocl.contexts[i],
					 CL_MEM_READ_WRITE,
					  (sizeof(int)*acc_distance_size),
					  NULL,&err);

    err = clSetKernelArg(ocl.kernel[i],
			 0,sizeof(cl_mem),ocl.acc_sequences+i);
    err = clSetKernelArg(ocl.kernel[i],
			 1,sizeof(cl_mem),ocl.acc_sequence_lengths+i);
    err = clSetKernelArg(ocl.kernel[i],
			 2,sizeof(cl_mem),ocl.acc_distances+i);
    err = clSetKernelArg(ocl.kernel[i],
			 3,sizeof(int),&seq_offset);
			       
  }
  free(gpu_sequence_dataset_buffer);

  ocl.num_devs =num_devs;
  ocl.devices = devices;
  ocl.platform = platform;
  ocl.contprop = contprop;
  ocl.programs = programs;
  

  //free(devices);
  //free(programs);
  
  return(ocl);
}

void opencl_destroy(opencl_stuff ocl) {
  free(ocl.cmdq);
  free(ocl.kernel);
  free(ocl.contexts);
  free(ocl.acc_sequences);
  free(ocl.acc_sequence_lengths);
  free(ocl.acc_distances);
  free(ocl.local_distances);
  free(ocl.programs);
  free(ocl.devices);
  free(ocl.contprop);
  free(ocl.platform);
}
  
#endif

static inline int compare(const void* a , const void* b) {
  int na = *(int*)a;
  int nb = *(int*)b;
  return (na-nb);
}

#if defined(_SCAN_SMITH_WATERMAN_GPU)
static inline neighbors region_query(int point, float epsilon, dataset ds,
				     opencl_stuff ocl) {

  int i,j,k;

  neighbors nb;

  cl_int err;
  
  size_t pointer_offset = 0;
  size_t acc_distance_size;
  
  nb.members=(int*)malloc(sizeof(int)*ds.n_values);

  for(i=0;i<ocl.num_devs;i++) {

    if (i == 0) {
      acc_distance_size = ds.n_values/ocl.num_devs + ds.n_values%ocl.num_devs;
    } else {
      acc_distance_size = ds.n_values/ocl.num_devs;
    }
    
    err = clSetKernelArg(ocl.kernel[i],4,sizeof(int),&point);
  
    err = clSetKernelArg(ocl.kernel[i],5,sizeof(int),&pointer_offset);
    
    err = clEnqueueNDRangeKernel(ocl.cmdq[i],ocl.kernel[i],1,
				 NULL,
				 &acc_distance_size,
				 NULL,
				 0,NULL,NULL);

    

    err = clEnqueueReadBuffer(ocl.cmdq[i], ocl.acc_distances[i],
			      CL_FALSE, 0,
			      acc_distance_size*sizeof(int),
			      ocl.local_distances+pointer_offset,
			      0,NULL,NULL);
    pointer_offset+=acc_distance_size;
  }
  
  for(i=0;i<ocl.num_devs;i++) {
    clFinish(ocl.cmdq[i]);
  }
  
  nb.n_members = 0;
  for(i =0;i<ds.n_values;i++) {
    if(ocl.local_distances[i] <= epsilon) {
      nb.members[nb.n_members++] = i;
    }
  }
  return(nb);
}

#elif defined(_SCAN_SMITH_WATERMAN)
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
  //qsort(nb.members,nb.n_members,sizeof(int),compare);
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
#if defined (_SCAN_SMITH_WATERMAN_GPU)
				  dataset ds, opencl_stuff ocl
#else				  
				  dataset ds
#endif
				  ) {

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
  memset(nb_unsorted.members,0,sizeof(int)*ds.n_values);
  
  merge = (int*)malloc(ds.n_values*sizeof(int));

  memcpy(nb_unsorted.members,nb.members,sizeof(int)*nb.n_members);
  nb_unsorted.n_members = nb.n_members;

  cl->members = (int*)malloc(sizeof(int)*ds.n_values);
  cl->members[cl->n_members] = point;
  cluster_member[point] = 1;
  cl->n_members++;

  for(i=0;i < nb_unsorted.n_members;i++) {
    if(!get_value_in_binary_array_at_index(visited,
					   (size_t)nb_unsorted.members[i])) {
      set_value_in_binary_array_at_index(visited,
					 (size_t)nb_unsorted.members[i]);
#if defined (_SCAN_SMITH_WATERMAN_GPU)
      nb_of_nb = region_query(nb_unsorted.members[i], epsilon, ds, ocl);
#else
      nb_of_nb = region_query(nb_unsorted.members[i], epsilon, ds);
#endif
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
      } else {
	free(nb_of_nb.members);
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
  free(merge);
}

split_set
#if defined (_SCAN_L1)
dbscan_L1
#elif defined (_SCAN_L2)
dbscan_L2
#elif defined (_SCAN_SMITH_WATERMAN)
dbscan_SW
#elif defined (_SCAN_SMITH_WATERMAN_GPU)
dbscan_SW_GPU
#endif
#if defined (_SCAN_SMITH_WATERMAN_GPU)
(dataset ds, float epsilon, int minpts, opencl_stuff ocl)
#else
(dataset ds, float epsilon, int minpts)
#endif
{

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
    if(!get_value_in_binary_array_at_index(visited,(size_t)i)) {
      set_value_in_binary_array_at_index(visited,(size_t)i);
#if defined (_SCAN_SMITH_WATERMAN_GPU)
      nb = region_query(i, epsilon, ds, ocl);
#else
      nb = region_query(i, epsilon, ds);
#endif
      if ( minpts < nb.n_members ) {
#if defined (_SCAN_SMITH_WATERMAN_GPU)
		expand_cluster(i,nb,clusters+(n_clusters),epsilon,minpts,
			       visited,cluster_member,ds, ocl);
#else
		expand_cluster(i,nb,clusters+(n_clusters),epsilon,minpts,
			       visited,cluster_member,ds);
#endif

		       
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

void adaptive_dbscan(
#if defined(_SCAN_SMITH_WATERMAN_GPU)
		     split_set (*dbscanner) (dataset,
					     float,
					     int,
					     opencl_stuff),
#else
		     split_set (*dbscanner) (dataset,
					     float,
					     int),
#endif
		     dataset ds,
		     float epsilon_start,
		     float epsilon_inc,
		     int minpts,
		     char* split_files_prefix
		     ) {
  
  int i,j,k;

  int initial_counter, count, eps_count;
  
  split_set* set_of_split_sets;
  cluster_connections** connections = NULL;
  cluster_connections* current_connection;

  cluster not_covered;

  float coverage;
  
  split_set new_split_set;

  char split_files[255];
  char buffer[20];

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  opencl_stuff ocl = opencl_initialization(ds);
#endif

  printf("Performing adaptive clustering with parameters: \n"
	 "minpts: %i\n"
	 "epsilon_start: %f\n"
	 "epsilon_inc: %f\n", minpts, epsilon_start, epsilon_inc);
  	 
  set_of_split_sets = (split_set*)malloc(sizeof(split_set));

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  set_of_split_sets[0] = dbscanner(ds, epsilon_start, minpts, ocl);
#else
  set_of_split_sets[0] = dbscanner(ds, epsilon_start, minpts);
#endif
  
  printf("Initial set obtained with %i clusters\n",
	  set_of_split_sets[0].n_clusters);

  initial_counter = 0;
  while( set_of_split_sets[0].n_clusters == 1 ) {
    if (initial_counter == 20 || epsilon_start == 0) {
      printf("Error did not find a sufficient" 
	     " starting position in 20 tries \n");
      _exit(1);
    }
    epsilon_start = epsilon_start/2;
    printf("Trying new starting point \n");

    free_split_set_and_associated_clusters(set_of_split_sets[0]);
#if defined(_SCAN_SMITH_WATERMAN_GPU)
    set_of_split_sets[0] = dbscanner(ds, epsilon_start, minpts, ocl);
#else
    set_of_split_sets[0] = dbscanner(ds, epsilon_start, minpts);
#endif
    initial_counter++;  
  }

  printf("Sarting with %i clusters\n", set_of_split_sets[0].n_clusters);
  
  not_covered = data_not_in_clusters(set_of_split_sets[0], ds);

  if(not_covered.n_members > 0) {
	free(not_covered.members);
  }

  printf("Coverage at initial point: %f\n",
	 (float)1.f-(float)not_covered.n_members/(float)ds.n_values);

  if(set_of_split_sets[0].n_clusters == 1) {
    printf("All clusters fusioned in one step, decrease epsilon increment\n");
    _exit(1);
  }

  count = 0;
  eps_count = 1;

  do{
#if defined(_SCAN_SMITH_WATERMAN_GPU)
    new_split_set =
      dbscanner(ds,epsilon_start+eps_count*epsilon_inc, minpts, ocl);
#else
    new_split_set = dbscanner(ds, epsilon_start+eps_count*epsilon_inc, minpts);
#endif
    if (new_split_set.n_clusters < set_of_split_sets[count].n_clusters) {
      
      not_covered = data_not_in_clusters(new_split_set, ds);
      
      printf("Layer %i has %i clusters\n",count+1, new_split_set.n_clusters);
      
      printf("  Coverage at layer %i: %f\n", count+1,
	     1.f-(float)not_covered.n_members/(float)ds.n_values);

      if(not_covered.n_members > 0) {
	free(not_covered.members);
      }

      printf("  Epsilon at layer %i: %f\n", count+1,
	     epsilon_start+eps_count*(float)epsilon_inc);
	          
      current_connection = generate_split_set_relation(set_of_split_sets[count],
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
    } else {
      free_split_set_and_associated_clusters(new_split_set);
    }
    eps_count++;
  }while(new_split_set.n_clusters != 1);

#if defined(_SCAN_SMITH_WATERMAN_GPU)
  opencl_destroy(ocl);
#endif
  
  printf("Connections found: \n");
  for(i=0;i<(count+1);i++) {
    
    sprintf(buffer,"%04d",i);
    memcpy(split_files,split_files_prefix,
	   strlen(split_files_prefix)+1);
    strcat(split_files,buffer);

    store_split_set(split_files, set_of_split_sets[i]);
  }

  for(i=count;i>0;i--) {
    printf("Layer %i:\n", i);
    for(j=0;j<set_of_split_sets[i].n_clusters;j++) {
      printf("  Cluster %i connected to: \n  ",j);
      for(k=0;k<connections[i-1][j].n_connections;k++) {
	printf("%i ", connections[i-1][j].connections[k]);
      }
      printf("\n");
    }
  }

  for(i=count;i>0;i--) {
    for(j=0;j<set_of_split_sets[i].n_clusters;j++) {
      free(connections[i-1][j].connections);
    }
    free(connections[i-1]);
  }
  free(connections);
    
  
  for(i=0;i<(count+1);i++) {
    free_split_set_and_associated_clusters(set_of_split_sets[i]);
  }
  free(set_of_split_sets);
}
  
  
