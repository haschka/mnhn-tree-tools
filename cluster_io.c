#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"dataset.h"
#include"cluster.h"
#include"binary_array.h"
#include"volumes.h"

double* array_deltas(double* array, int array_length) {

  int i;
  
  double* diffs = (double*)malloc(sizeof(double)*(array_length-1));

  for(i=0;i<array_length-1;i++) {
    diffs[i] =  array[i] - array[i+1];
  }
  return(diffs);
}

double* get_epsilon_dist_from_adaptive_clustering_output(FILE* f,
							 int n_layers) {
 
  int i;
  char* line = NULL;
  size_t line_size;
                         
  char inital_check[] = "epsilon_start:"; /* 13 chars */
  char layer_check[] =  "  Epsilon at layer"; /* 17 chars */

  double * epsilons = (double*)malloc(sizeof(double)*n_layers);
  int epsilon_counter = 0;
  
  rewind(f);
  
  while ( -1 != getline(&line, &line_size, f)) {

    if (line_size > 10) {

      if (epsilon_counter == 0) {
      
	for(i=0;i<14;i++) {
	  if (inital_check[i] != line[i]) {
	    goto check_layer;
	  }
	}
	sscanf(line,"%*s %lf", epsilons);
	epsilon_counter = 1;
      }
      
    check_layer:
      
      for(i=0;i<17;i++) {
	if(layer_check[i] != line[i]) {
	  goto checks_finished;
	}
      }
      
      if(epsilon_counter < 1) {
	printf("Malformed output file epsilon_start not found before \n "
	       "  Epsilon at layer: \n");
	_exit(1);
      }

      if(epsilon_counter == n_layers) {
	goto all_read;
      }
      
      sscanf(line,"  %*s %*s %*s %*s %lf", epsilons+epsilon_counter);
      epsilon_counter++;
      
    checks_finished:
      i=0;
      
    }
  }
  
 all_read:
  if (epsilon_counter < n_layers) {
    printf("Not enough layers in output file! \n");
    printf("Counted %i layers, requested %i layers! \n",
	   epsilon_counter,
	   n_layers);
    _exit(1);
  }
  free(line);
  return(epsilons);
}

double* densities_from_epsilons(int type, int dimensions,
				double* epsilons, int n_epsilons,
				int min_points) {

  int i;

  double(*volume_function)(int,double);

  double* densities = (double*)malloc(sizeof(double)*n_epsilons);
  
  switch(type) {
  case 1:
    volume_function = &vol_hypercube;
    break;
  case 2:
    volume_function = &vol_hypersphere;
    break;
  default:
    printf("Wrong parameter to densities_from_epsilon (min_points)\n");
    _exit(1);
  }
  
  for (i=0;i<n_epsilons;i++) {
    densities[i] = min_points/(volume_function(dimensions,epsilons[i]));
  }

  return(densities);
}							      

tree_node* generate_tree(int n_layers, cluster_connections** c,
			 split_set *sets, double *tree_lengths) {
  
  int i, j, k;

  tree_node* current_node;
  
  tree_node** nodes = (tree_node**)malloc(sizeof(tree_node*)*(n_layers));
  
  for(i=0;i<n_layers;i++) {
    nodes[i] = (tree_node*)malloc(sizeof(tree_node)
				  *sets[i].n_clusters);
    for(j=0;j<sets[i].n_clusters;j++) {
      nodes[i][j].id = NULL;
      if(tree_lengths != NULL) {
	if(i==(n_layers-1)) {
	  nodes[i][j].length = 0;
	} else {
	  nodes[i][j].length = tree_lengths[i];
	}
      } else {
	if(i==(n_layers-1)) {
	  nodes[i][j].length = 0;
	} else {
	  nodes[i][j].length = 1;
	}
      }	  
      nodes[i][j].n_members = 0;
      nodes[i][j].child = NULL;
      nodes[i][j].neighbor = NULL;
      nodes[i][j].parent = NULL;
    }
  }
  
  for(i=(n_layers-1);i>0;i--) {
    for(j=0;j<sets[i].n_clusters;j++) {
      for(k=0;k<c[i-1][j].n_connections;k++) {
	current_node = nodes[i-1]+c[i-1][j].connections[k];
	current_node->parent = nodes[i]+j;
	if(k>0) {
	  current_node->neighbor =
	    nodes[i-1]+c[i-1][j].connections[k-1];
	} else {
	  current_node->neighbor = NULL;
	}
	current_node->id = (char*)malloc(sizeof(char)*20);
	current_node->n_members =
	  sets[i-1].clusters[c[i-1][j].connections[k]].n_members;
	sprintf(current_node->id,"L%iC%iN%i",
		i-1,
		c[i-1][j].connections[k],
		current_node->n_members);

      }
      if(c[i-1][j].n_connections > 0) {
	nodes[i][j].child =
	  nodes[i-1]+c[i-1][j].connections[c[i-1][j].n_connections-1];
      }
    }
  }

  nodes[n_layers-1]->id = (char*)malloc(sizeof(char)*20);
  sprintf(nodes[n_layers-1]->id,"L%iC%i",n_layers-1,0);
  
  return(nodes[n_layers-1]);
}

split_set filtered_split_set_by_min_size(split_set s_in, int min_size) {

  int i, count, current_n_members;
  
  split_set s_out;

  s_out.clusters = (cluster*)malloc(sizeof(cluster)*s_in.n_clusters);
  
  count = 0;
  for(i=0;i<s_in.n_clusters;i++) {
    current_n_members = s_in.clusters[i].n_members;
    if ( current_n_members >= min_size) {
      s_out.clusters[count].n_members = current_n_members;
      s_out.clusters[count].members =
	(int*)malloc(sizeof(int)*current_n_members);
      memcpy(s_out.clusters[count].members,
	     s_in.clusters[i].members,
	     sizeof(int)*current_n_members);
      count++;
    }
  }
  s_out.n_clusters = count;
  s_out.clusters = (cluster*)realloc(s_out.clusters,sizeof(cluster)*count);
  return(s_out);
}

void print_tree_worker(FILE*f, tree_node* root)
{
  tree_node* n;
  if(root->child!=NULL) {
    fprintf(f,"(");
    for(n=root->child;n!=NULL;n=n->neighbor) {
      if(n!=root->child) fprintf(f,",");
      print_tree_worker(f,n);
    }
    fprintf(f,")");
  }
  fprintf(f,"%s",root->id);
  if(root->parent!=NULL) printf(":%f",root->length);
}

void print_tree(FILE*f, tree_node* root) {
  print_tree_worker(f,root);
  fprintf(f,";");
}

double pureness_helper(int layer, int this_cluster,
		       long long int* inpures,
		       split_set *s, split_set target) {

  int i;
  
  cluster current_cluster = s[layer].clusters[this_cluster];
  cluster intersection;

  double inpureness = 0;

  int c;
  
  for(i = 0; i < target.n_clusters; i++) {
    intersection =
      intersection_of_clusters(current_cluster, target.clusters[i]);
    if(current_cluster.n_members != intersection.n_members) {
      inpures[layer]++;
      if ( intersection.n_members > current_cluster.n_members/2 ) {
	c = current_cluster.n_members - intersection.n_members;
      } else {
	c = intersection.n_members;
      }
      inpureness+=(double)c/((double)current_cluster.n_members-(double)c);
    }
  }
  return(inpureness);
}
  
void pureness_worker(double* pureness, long long int * inpures,
		     tree_node* root, split_set *s,
		     split_set target) {
  tree_node* n;
  int layer, this_cluster;
  if(root->child!=NULL) {
    for(n=root->child;n!=NULL;n=n->neighbor) {
      pureness_worker(pureness,inpures,n,s,target);
    }
  }
  sscanf(root->id,"L%iC%i", &layer, &this_cluster);
  pureness[layer] += pureness_helper(layer, this_cluster, inpures, s, target);
}
  
double* pureness_from_tree(int n_layers, tree_node* root,
			   split_set *s, split_set target) {

  double *pureness = (double*)malloc(sizeof(double)*n_layers);
  long long int* inpures = (long long*)malloc(sizeof(long long int)*n_layers);

  int i;
  
  memset(pureness,0,sizeof(double)*n_layers);
  memset(inpures,0,sizeof(long long int)*n_layers);
  
  pureness_worker(pureness, inpures, root, s, target);
  
  for(i = 0 ; i < n_layers ; i++) {
    pureness[i] =
      pureness[i]/(double)inpures[i];
  }
  free(inpures);
  return(pureness);
}

double* clusters_in_layer_vs_target_clusters(int n_layers, split_set *s,
					     split_set target) {

  int i;
  double *cluster_factors = (double*)malloc(sizeof(double)*n_layers);

  for(i = 0 ; i < n_layers ; i++) {
    cluster_factors[i] = ((double)(abs(s[i].n_clusters-target.n_clusters)))
      / (double)(target.n_clusters);
  }
  return(cluster_factors);
}

void create_single_cluster_file_with_values(char* filename, cluster cl,
					    dataset ds) {
  int j, k;
  FILE* f = fopen(filename, "w");

  for (j = 0; j< cl.n_members; j++) {
    for ( k = 0; k < ds.n_dimensions - 1; k++) {
      fprintf(f, "%f\t", ds.values[k][cl.members[j]]);
    }
    fprintf(f, "%f\n",
	    ds.values[ds.n_dimensions-1][cl.members[j]]);
  }
  fclose(f);
}

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

    create_single_cluster_file_with_values(file_name, s.clusters[i], ds);
    
  }
}

static inline int comp_int(const void *a, const void* b) {
  int* a_i = (int*)a;
  int* b_i = (int*)b;

  return(a_i[0]-b_i[0]);
}

static inline int binary_search(int* array, size_t array_size, int target) {

  int left = 0;
  int right = array_size - 1;
  int pos;
  
  while( left <= right ) {
    pos = (left+right)/2;
    if (array[pos] < target) {
      left = pos + 1;
    } else if (array[pos] > target) {
      right = pos - 1;
    } else {
      return pos;
    }
  }
  return -1;
}

cluster_connections*
generate_split_set_relation(split_set ancient, split_set new) {

  int i,j, counter;
  
  cluster_connections* connections_ancient_to_new =
    (cluster_connections*)malloc(sizeof(cluster_connections)*new.n_clusters);
  
  cluster intersection;
  
  for(i=0;i<new.n_clusters;i++) {

    connections_ancient_to_new[i].connections =
      (int*)malloc(sizeof(int)*ancient.n_clusters);
    counter = 0;
    for(j=0;j<ancient.n_clusters;j++) {

      intersection = intersection_of_clusters(ancient.clusters[j],
					      new.clusters[i]);

      if ( (double)intersection.n_members >
	   (double)ancient.clusters[j].n_members*(double)0.8 ) {

	connections_ancient_to_new[i].connections[counter] = j;
	counter++;

      }
      if (intersection.n_members > 0) {
	free(intersection.members);
      }
    }
    connections_ancient_to_new[i].connections =
      (int*)realloc(connections_ancient_to_new[i].connections,
		    sizeof(int)*counter);
    connections_ancient_to_new[i].n_connections = counter;
  }
  return(connections_ancient_to_new);
}

cluster intersection_of_clusters(cluster a, cluster b) {

  cluster big;
  cluster small;

  int i, j, count;
  int search_result;
  
  cluster intersection;
  
  if ( a.n_members >= b.n_members ) {
    big = a;
    small = b;
  } else {
    big = b;
    small = a;
  }

  intersection.members= (int*)malloc(sizeof(int)*small.n_members);

  count = 0;
  if( big.n_members > 20 ) {
  
    qsort(big.members, big.n_members, sizeof(int), comp_int);

    for(i = 0 ; i < small.n_members; i++) {
      
      search_result = binary_search(big.members,
				    big.n_members,
				    small.members[i]);

      if (search_result != -1) {
	intersection.members[count] = small.members[i];
	count++;
      }
    }
    
  } else {

    for(i = 0 ; i < small.n_members; i++) {
      for(j = 0; j < big.n_members;j++) {
	if (small.members[i] == big.members[j]) {
	  intersection.members[count] = small.members[i];
	  count++;
	}
      }
    }
    
  }
  intersection.n_members = count;
  intersection.id = -1;
  intersection.members = (int*)realloc(intersection.members, sizeof(int)*count);

  return(intersection);
}

void print_cluster_matrix_view_annotation(FILE* f, dataset ds) {

  int i;

  for(i=0;i<ds.n_values-1;i++) {
    fprintf(f,"%i\t",i);
  }
  fprintf(f,"%i\n",ds.n_values-1);
}

void print_cluster_matrix_view(FILE*f, split_set s, dataset ds) {

  int i,j;

  int* matrix_line = (int*)malloc(sizeof(int)*ds.n_values);

  for(i=0;i<ds.n_values;i++) {
    matrix_line[i] = -1;
  }

  for( i = 0; i< s.n_clusters; i++) {
    for (j = 0; j< s.clusters[i].n_members; j++) {
      matrix_line[s.clusters[i].members[j]] = i;
    }
  }

  for( i = 0; i<ds.n_values-1;i++) {
    fprintf(f,"%i\t", matrix_line[i]);
  }
  fprintf(f,"%i\n", matrix_line[ds.n_values-1]);
  free(matrix_line);
}

cluster cluster_from_sequence_in_dataset(dataset ds,
					 char* seq,
					 size_t seq_len) {

  int i,j;

  cluster cl;

  cl.members = (int*)malloc(sizeof(int)*ds.n_values);
  cl.n_members = 0;

  i = 0;
  while(i<ds.n_values) {
    if (ds.sequence_lengths[i] == seq_len) {
      for(j=0;j<seq_len;j++) {
	if(ds.sequences[i][j] != seq[j]) {
	  goto sequence_false;
	}
      }
      cl.members[cl.n_members++] = i;
    }
  sequence_false:
    i++;
  }
  cl.members = (int*)realloc(cl.members,sizeof(int)*cl.n_members);
  return(cl);
}

cluster data_not_in_clusters(split_set s, dataset ds) {

  int i,j;
  
  cluster out;

  char* out_b = alloc_and_set_zero_binary_array(ds.n_values);
  size_t count; 

  out.members = (int*)malloc(sizeof(int)*ds.n_values);
  
  for( i = 0; i< s.n_clusters; i++) {
    for (j = 0; j< s.clusters[i].n_members; j++) {
      set_value_in_binary_array_at_index( out_b, s.clusters[i].members[j]);
    }
  }

  count = 0;
  for( i = 0; i< ds.n_values; i++) {
    if(!(get_value_in_binary_array_at_index(out_b,i))){
      out.members[count] = i;
      count++;
    }
  }
  out.id = -1;
  out.members = (int*)realloc(out.members,sizeof(int)*count);
  out.n_members = count;
  free(out_b);
  return(out);
}

split_set read_split_set(char* filename) {

  char s_set[] = "SPLIT_SET";
  char buffer[10];
  
  int i,j;

  int fd = open(filename, O_RDONLY);

  size_t check = 0;
  size_t count;

  split_set s;

  size_t si = sizeof(int);
  
  check += read(fd,buffer,9);
  buffer[9] = 0;
  if (0 != strcmp(buffer,s_set)) {
    printf("Error %s is not a split_set\n", filename);
    _exit(1);
  }
  check += read(fd, &s.n_clusters, sizeof(int));
  count = 9+si;
  s.clusters = (cluster*)malloc(sizeof(cluster)*s.n_clusters);
  for(i = 0; i < s.n_clusters; i++) {
    check += read(fd, &s.clusters[i].n_members, sizeof(int));
    check += read(fd, &s.clusters[i].id, sizeof(int));
    count += 2*si;
    s.clusters[i].members = (int*)malloc(sizeof(int)*s.clusters[i].n_members);
    for(j= 0; j < s.clusters[i].n_members; j++) {
      check += read(fd, s.clusters[i].members+j, sizeof(int));
      count += si;
    }
  }
  if (check != count) {
    printf("There was an error reading " 
	   "the binary clusterfile %s", filename);
  }
  close(fd);
  return(s);
}   
  
void store_split_set(char* filename, split_set s) {

  char s_set[] = "SPLIT_SET";

  int i,j;
  
  int fd = open(filename, O_WRONLY | O_CREAT,
		S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);

  size_t check = 0;
  size_t count;
  
  size_t si = sizeof(int);
  
  check += write(fd, s_set, 9);
  check += write(fd, &s.n_clusters, sizeof(int));
  count = 9+si;
  for(i = 0; i < s.n_clusters; i++) {
    check += write(fd, &s.clusters[i].n_members, sizeof(int));
    check += write(fd, &s.clusters[i].id, sizeof(int));
    count += 2*si;
    for(j= 0; j < s.clusters[i].n_members; j++) {
      check += write(fd, s.clusters[i].members+j, sizeof(int));
      count += si;
    }
  }
  if (check != count) {
    printf("There was an error in creating " 
	   "the binary clusterfile %s", filename);
  }
  close(fd);
}
  

void create_single_cluster_file(char* filename, cluster cl, dataset ds) {

  int j, k;
  FILE* f = fopen(filename, "w");

  for (j = 0; j< cl.n_members; j++) {
    fprintf(f,">sequence_%i\n",cl.members[j]);
    for ( k = 0;
	  k < (strlen(ds.sequences[cl.members[j]])-1);
	  k++) {
      if( k != 0 && k%50 == 0 ) {
	fprintf(f, "\n");
	fputc(ds.sequences[cl.members[j]][k],f);
      } else {
	fputc(ds.sequences[cl.members[j]][k],f);
      }
    }
    if ( (strlen(ds.sequences[cl.members[j]])-1) != 0 &&
	 (strlen(ds.sequences[cl.members[j]])-1) %50 == 0 ) {
      fprintf(f, "\n");
      fputc(ds.sequences[cl.members[j]][k],f);
    } else {
      fputc(ds.sequences[cl.members[j]][k],f);
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

void create_cluster_files(char* prefix, split_set s, dataset ds) {

  int i;

  char file_name[256];

  char i_buffer[6];
  
  for ( i = 0; i< s.n_clusters; i++) {

    file_name[0] = 0;
    strcat(file_name, prefix);
    strcat(file_name, "-");
    sprintf(i_buffer,"%03d", i);
    strcat(file_name, i_buffer);

    create_single_cluster_file(file_name,s.clusters[i],ds);
    
  }
}

void free_split_set_and_associated_clusters(split_set s) {

  int i;
  for(i=0;i<s.n_clusters;i++) {
    free(s.clusters[i].members);
  }
  free(s.clusters);
}
    
