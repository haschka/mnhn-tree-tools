#include<stdio.h>


int main(int argc, char* argv) {
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);
  
  if ( NULL == (projection_f = fopen(argv[2], "r"))) file_error(argv[2]);
  load_projections_from_file_into_dataset(projection_f,dimensions,&ds);
  fclose(projection_f);
