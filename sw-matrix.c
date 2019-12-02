#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"dataset.h"
#include"smith-waterman.h"

int main(int argc, char** argv) {

  FILE* fasta_one;
  FILE* fasta_two;

  dataset ds_one;
  dataset ds_two;

  size_t i,j;

  int* matrix;

  double sum;
  double mean;
  double sigma;
  
  fasta_one = fopen(argv[1],"r");

  ds_one = dataset_from_fasta(fasta_one);
  fclose(fasta_one);

  fasta_two = fopen(argv[2],"r");

  ds_two = dataset_from_fasta(fasta_two);

  matrix = (int*)malloc(sizeof(int)*ds_one.n_values*ds_two.n_values);
  
  for(i=0;i<ds_one.n_values;i++) {
    for(j=0;j<ds_two.n_values;j++) {

      matrix[i*ds_two.n_values+j] =
	score(ds_one.sequences[i],ds_two.sequences[j]);

      printf("%lu %lu %u \n", i, j, matrix[i*ds_two.n_values+j]);
    }
  }

  sum = 0.;
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    sum+=(double)matrix[i];
  }
  printf("sum = %lf, ",sum);
  mean = sum/((double)ds_one.n_values*(double)ds_two.n_values);

  sigma = 0.;
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    sigma+=((double)matrix[i]-(double)mean)*((double)matrix[i]-(double)mean);
  }
  sigma = sqrt(sigma)/(ds_one.n_values+ds_two.n_values);

  printf("mean = %lf, sigma = %lf \n", mean, sigma);
 
}

  
  
