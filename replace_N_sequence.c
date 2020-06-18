#include<stdio.h>
#include<unistd.h>

#include"dataset.h"
#include"kmers.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

void print_arguments() {

  printf("Arguments are: \n"
	 "   [file-in ] FASTA file with a sequence containing Ns \n"
	 "   [int     ] n_th sequence in FASTA input file to be treated \n"
	 "   [file-out] FASTA file that will contain possible replacements \n"
	 );
}

int main(int argc, char** argv) {

  FILE* f_in;
  FILE* f_out;
  
  dataset ds_in;
  dataset ds_out;

  unsigned long idx;

  if(argc < 3) {
    print_arguments();
    return(1);
  }
  
  sscanf(argv[2], "%lu", &idx);

  if ( NULL == (f_in = fopen(argv[1], "r"))) file_error(argv[1]);
  ds_in = dataset_from_fasta(f_in);
  fclose(f_in);

  ds_out = generate_n_permutations_for_sequence(ds_in, (size_t)(idx-1));
  free_sequences_from_dataset(ds_in);
  
  if ( NULL == (f_out = fopen(argv[3], "w"))) file_error(argv[3]);
  dataset_to_fasta(f_out,ds_out);
  fclose(f_out);
  free_sequences_from_dataset(ds_out);

  return(0);
}
  
  
  
  

  
    
