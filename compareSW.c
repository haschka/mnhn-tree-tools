#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"dataset.h"
#include"comparison.h"

int main(int argc, char** argv) {

  FILE* fasta_one;
  FILE* fasta_two;

  dataset ds_one;
  dataset ds_two;

  size_t i,j;

  unsigned long* matrix;

  double mean;
  double sigma;

  size_t n_threads;

  if (argc < 3) {
    printf("Arguments are [fasta1] [fasta2] [n_threads]\n");
    return(1);
  }
  
  fasta_one = fopen(argv[1],"r");

  ds_one = dataset_from_fasta(fasta_one);
  fclose(fasta_one);

  fasta_two = fopen(argv[2],"r");

  ds_two = dataset_from_fasta(fasta_two);
  fclose(fasta_two);

  sscanf(argv[3],"%lu",&n_threads);
  
  mean = mean_from_smith_waterman_datasets(ds_one,ds_two,n_threads);
  sigma = sigma_from_smith_waterman_datasets( mean,
					      ds_one, ds_two,
					      n_threads);
  printf("mean = %lf, sigma = %lf \n", mean, sigma);
 
}

  
  
