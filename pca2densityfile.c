#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include"dataset.h"
#include"density.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {

  dataset ds;

  size_t input_dimensions;
  int output_dimensions = 2;
  
  FILE* fasta_f;
  FILE* projection_f;

  density_map lmap, map;

  float interval_length;

  if(argc < 4) {
    printf("Arguments are: \n"
	   " [file] Sequences in FASTA \n"
	   " [file] Corresponding PCA from kmers \n"
	   " input-dimensions: dimensions in projection file \n"
	   " interval-length for density pixel/voxel integration \n"
	   );
    return(1);
  }
  
  sscanf(argv[3], "%lu", &input_dimensions);
  sscanf(argv[4], "%f", &interval_length);
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (projection_f = fopen(argv[2], "r"))) file_error(argv[2]);
  load_projections_from_file_into_dataset(projection_f,input_dimensions,&ds);
  fclose(projection_f);

  lmap = generate_density_map_from_dataset(ds, output_dimensions,
					  interval_length);

  printf("   shift x:  %4.8f\n", lmap.shift[0]);
  printf("   shift y:  %4.8f\n", lmap.shift[1]);
  printf("n points x:  %5lu\n", lmap.n_values_per_dimension[0]);
  printf("n points y:  %5lu\n", lmap.n_values_per_dimension[1]);
  print_2d_density_to_file(stdout, lmap);

  free(lmap.intensities);
  return(0);
}
  
