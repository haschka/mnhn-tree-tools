#include<stdio.h>
#include<unistd.h>
#include"dataset.h"
#include"comparison.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

int main(int argc, char** argv) {


  FILE* fasta_f;
  FILE* projection_f;
  
  int dimensions;
  dataset ds;

  if(argc < 3) {
    printf("Arguments are: \n");
    printf("  [FILE-in] fasta file to print distance in different norms from");
    printf("  [FILE-in] PCA projections file\n");
    printf("  (int) number of projections in projections file\n");
    _exit(1);
  }
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  sscanf(argv[3], "%i", &dimensions);
  
  if ( NULL == ( projection_f = fopen(argv[2], "r"))) file_error(argv[2]);
  load_projections_from_file_into_dataset(projection_f,dimensions,&ds);
  fclose(projection_f);

  print_SW_PCA_L_comparison(stdout, ds);
}
