#include<stdio.h>

#include"dataset.h"
#include"comparison.h"

int main(int argc, char** argv) {

  FILE* f_one;
  FILE* f_two;

  dataset ds_one, ds_two;
  
  f_one = fopen(argv[1], "r");

  ds_one = dataset_from_fasta(f_one);
  fclose(f_one);

  f_two = fopen(argv[2], "r");
  
  ds_two = dataset_from_fasta(f_two);
  fclose(f_two);

  printf("%lf\n", silhouette_from_smith_waterman_datasets(ds_one,
							  ds_two));
}
	 

    
  
