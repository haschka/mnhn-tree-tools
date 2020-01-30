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
  int output_dimensions;
  
  FILE* fasta_f;
  FILE* projection_f;

  density_map map;

  float interval_length;

  if(argc < 6) {
    printf("Arguments are: \n"
	   " [file] Sequences in FASTA \n"
	   " [file] Corresponding PCA from kmers \n"
	   " input-dimensions: dimensions in projection file \n"
	   " output-dimensions: 2 or 3 for density map \n"
	   " interval-length for density pixel/voxel integration \n"
	   " [file/path-prefix] filename or path/prefix for density png files\n"
	   );
    return(1);
  }
  
  sscanf(argv[3], "%lu", &input_dimensions);
  sscanf(argv[4], "%i", &output_dimensions);
  sscanf(argv[5], "%f", &interval_length);
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (projection_f = fopen(argv[2], "r"))) file_error(argv[2]);
  load_projections_from_file_into_dataset(projection_f,input_dimensions,&ds);
  fclose(projection_f);

  map = generate_density_map_from_dataset(ds, output_dimensions,
					  interval_length);

  switch(output_dimensions) {
  case 2:
    save_2d_density_to_png(argv[6], map);
  case 3:
    save_3d_density_to_pngs(argv[6], map);
  default:
    return(1);
  }
  free(map.intensities);
  return(0);
}
  
