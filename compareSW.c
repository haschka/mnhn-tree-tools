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
  
  fasta_one = fopen(argv[1],"r");

  ds_one = dataset_from_fasta(fasta_one);
  fclose(fasta_one);

  fasta_two = fopen(argv[2],"r");

  ds_two = dataset_from_fasta(fasta_two);
  fclose(fasta_two);
  
  matrix = smith_waterman_distances_matrix(ds_one,ds_two);
  mean = mean_from_smith_waterman_distance_matrix(matrix,ds_one,ds_two);
  sigma = sigma_from_smith_waterman_distance_matrix(matrix, mean,
						    ds_one, ds_two);
  printf("mean = %lf, sigma = %lf \n", mean, sigma);
 
}

  
  
