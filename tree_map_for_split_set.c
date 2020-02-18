#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#include"dataset.h"
#include"cluster.h"

typedef struct {
  int cluster_index;
  int r;
  int g;
  int b;
} cluster_color;

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

void print_arguments() {
  printf("Arguments are: \n"
	 "  [split_set] file containing clusters to be mapped onto tree\n"
	 "  [fasta] database of sequences the split sets are refering to\n"
	 "  (int) number of clusters to color from splitset \n"
	 "  (int) 0 or 1 use gradient for number of sequences from cluster \n"
	 "        if 1 the colors are given from: white,\n"
	 "                                         1 sequence \n"
	 "                                        to the color selected,\n"
	 "                                         all sequences \n"
	 
	 "  [(int),string ... (int),string] number of cluster to color in \n"
	 "                                  in the color given by the "
	 "string \n"
	 "                                  #1122FF\n"
	 "                                  number of tuples has to \n"
	 "                                  correspond to the number of \n"
	 "                                  clusters ( prev. argument ) \n"
	 "  [split_sets...] cluster files corresponding to newick tree\n");
}

cluster_color map_input_string_to_cluster_color(char* in, int max_clusters) {

  int i;
  int cluster_index;
  int pos_virgule = 0;
  char colorbuffer[256];
  char indexbuffer[256];

  char color_conversion_buffer[3];
  
  cluster_color ret_val;
  
  size_t in_length = strnlen(in, 256);
  size_t color_length;
  
  for (i = 0 ;i<in_length;i++) {
    if(in[i] == ',') {
      pos_virgule = i;
    }
  }
  if (pos_virgule == 0) {
    print_arguments();
  }
  for(i=0 ;i<pos_virgule;i++) {
    indexbuffer[i] = in[i];
  }

  if(pos_virgule == 0) {
    printf("Input Malformatted! \n");
    print_arguments();
    _exit(1);
  }
  
  indexbuffer[pos_virgule] = 0;
  for(i=pos_virgule+1;i<in_length;i++) {
    colorbuffer[i-(pos_virgule+1)]= in[i];
  }
  colorbuffer[in_length-(pos_virgule+1)] = 0;

  sscanf(indexbuffer, "%i",&cluster_index);

  if(cluster_index < 1 || cluster_index > max_clusters) {
    printf("Input malformatted cluster_index out of range %i not in [1-%i] \n",
	   cluster_index, max_clusters);
    print_arguments();
    _exit(1);
  }
  color_length = strnlen(colorbuffer,255);
  if(color_length != 7) {
    printf("Input malformatted color not of correct length \n!");
    print_arguments();
    _exit(1);
  }
  
  ret_val.cluster_index = cluster_index - 1;

  /* red */
  for(i=1;i<3;i++) color_conversion_buffer[i-1] = colorbuffer[i];
  color_conversion_buffer[3] = 0;
  ret_val.r = strtol(color_conversion_buffer,NULL, 16);

  /* green */
  for(i=3;i<5;i++) color_conversion_buffer[i-3] = colorbuffer[i];
  color_conversion_buffer[3] = 0;
  ret_val.g = strtol(color_conversion_buffer,NULL, 16);

  /* blue */
  for(i=5;i<7;i++) color_conversion_buffer[i-5] = colorbuffer[i];
  color_conversion_buffer[3] = 0;
  ret_val.b = strtol(color_conversion_buffer,NULL, 16);

  return(ret_val);
}
  
static inline int is_cluster_selected(cluster_color* cl_c, int current,
				      int n_input_clusters) {
  int i;
  int is_true = 0;
  for(i=0;i<n_input_clusters;i++) {
    if( cl_c[i].cluster_index == current ) {
      is_true = 1;
    }
  }
  return(is_true);
}

cluster_color get_color_from_cluster_index(cluster_color* cl_c,
					   int current,
					   int n_input_clusters) {
  int i;
  for(i=0;i<n_input_clusters;i++) {
    if( cl_c[i].cluster_index == current) {
      return(cl_c[i]);
    }
  }
}
  

int main(int argc, char** argv) {

  int i, j, k,l;
  
  FILE* fasta_s;
  FILE* fasta_f;

  dataset ds_search;
  dataset ds;
  split_set s;
  split_set input_split_set;
  
  cluster cl_search;
  cluster intersection;

  int n_split_sets;

  int n_input_clusters;
  int gradient_bool;

  float fraction_in_cluster;
  
  cluster_color* cluster_colors;

  cluster_color current_color;

  char out_color[8];
  
  if (argc < 5) {
    print_arguments();
    return(1);
  }
  
  input_split_set = read_split_set(argv[1]);

  if ( NULL == (fasta_f = fopen(argv[2], "r"))) file_error(argv[2]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  sscanf(argv[3],"%i",&n_input_clusters);
  
  if(n_input_clusters < 1 || n_input_clusters > input_split_set.n_clusters) {
    printf("Number of cluster incoherent: split_set contains %i clusters!\n",
	   input_split_set.n_clusters);
    print_arguments();
  }

  sscanf(argv[4],"%i",&gradient_bool);
  
  if(gradient_bool != 1 && gradient_bool != 0) {
    printf("Malformed input for gradient %i is not 1 or 0\n",gradient_bool);
    print_arguments();
    return(1);
  }

  cluster_colors =
    (cluster_color*)malloc(sizeof(cluster_color)*n_input_clusters);

  for(i=5;i<(5+n_input_clusters);i++) {
    cluster_colors[i-5] =
      map_input_string_to_cluster_color(argv[i],
					input_split_set.n_clusters);
  }

  n_split_sets = argc - (5+n_input_clusters);
  
  for(i=0;i<n_split_sets;i++) {
    s = read_split_set(argv[i+5+n_input_clusters]);
    for(j= 0;j<s.n_clusters;j++) {
      for(k= 0; k<input_split_set.n_clusters;k++) {
	if(is_cluster_selected(cluster_colors, k, n_input_clusters)) {
	  intersection = intersection_of_clusters(s.clusters[j],
						  input_split_set.clusters[k]);
	  if(intersection.n_members > 0) {
	    current_color = get_color_from_cluster_index(cluster_colors,
							 k,
							 n_input_clusters);
	    if(gradient_bool) {
	      fraction_in_cluster =
		(float)intersection.n_members
		/(float)input_split_set.clusters[k].n_members;
	      current_color = get_color_from_cluster_index(cluster_colors,
							   k,
							   n_input_clusters);
	      current_color.r =
		(int)((float)current_color.r*(float)fraction_in_cluster);
	      current_color.g =
		(int)((float)current_color.g*(float)fraction_in_cluster);
	      current_color.b =
		(int)((float)current_color.b*(float)fraction_in_cluster);
	    }
	    sprintf(out_color,"#%02x%02x%02x",
		    current_color.r,current_color.g,current_color.b);

	    for(l=1;l<8;l++) {
	      switch(out_color[l]) {
	      case 'a': out_color[l] = 'A'; break;
	      case 'b': out_color[l] = 'B'; break;
	      case 'c': out_color[l] = 'C'; break;
	      case 'd': out_color[l] = 'D'; break;
	      case 'e': out_color[l] = 'E'; break;
	      case 'f': out_color[l] = 'F'; break;
	      default: break;
	      }
	    }
	    
	    
	    printf("stroke:%s I ",out_color);
	    printf("L%iC%iN%i \n",i,j,s.clusters[j].n_members);
	  }
	}
      }
    }
    free_split_set_and_associated_clusters(s);
  }
  free(cluster_colors);
  free_sequences_from_dataset(ds);
  free_split_set_and_associated_clusters(input_split_set);
  return(0);
}
    
  
