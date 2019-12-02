#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<limits.h>
#include"dataset.h"
#include"smith-waterman.h"

unsigned long* smith_waterman_distances_matrix(dataset ds_one,
						dataset ds_two) {

  unsigned long* matrix =
    (unsigned long*)malloc(sizeof(unsigned long)
			   *ds_one.n_values*ds_two.n_values);

  size_t i,j;
  
  for(i=0;i<ds_one.n_values;i++) {
    for(j=0;j<ds_two.n_values;j++) {

      matrix[i*ds_two.n_values+j] =
	score(ds_one.sequences[i],ds_two.sequences[j]);
    }
  }

  return(matrix);
}

double mean_from_smith_waterman_datasets(dataset ds_one, dataset ds_two) {

  size_t i,j;

  double sum = 0.;
  
  for(i=0;i<ds_one.n_values;i++) {
    for(j=0;j<ds_two.n_values;j++) {

      sum+=
	score(ds_one.sequences[i],ds_two.sequences[j]);
    }
  }
  return(sum/((double)ds_one.n_values*(double)ds_two.n_values));
}

double sigma_from_smith_waterman_datasets(double mean,
					  dataset ds_one,
					  dataset ds_two) {

  size_t i,j;
  double sigma = 0.;
  double v;
  
  for(i=0;i<ds_one.n_values;i++) {
    for(j=0;j<ds_two.n_values;j++) {
      v = score(ds_one.sequences[i],ds_two.sequences[j]);
      sigma+=((double)v-(double)mean)*((double)v-(double)mean);
    }
  }
  return(sqrt(sigma)/(ds_one.n_values+ds_two.n_values));
}

double silhouette_from_smith_waterman_datasets(dataset ds_one,
					       dataset ds_two) {

  double internal;
  double external;

  double a,b,s, s_mean;
  
  size_t i,j;

  s_mean = 0;
  for(i=0;i<ds_one.n_values;i++) {
    internal = 0;
    external = 0;
    for(j=0;j<ds_one.n_values;j++) {
      internal += (double)score(ds_one.sequences[i],ds_one.sequences[j]);
    }
    for(j=0;j<ds_two.n_values;j++) {
      external += (double)score(ds_one.sequences[i],ds_two.sequences[j]);
    }
    
    a=(internal/((double)ds_one.n_values-1.));
    b=(external/((double)ds_two.n_values));

    if(a<b) {
      s = 1-a/b;
    } else if(a==b) {
      s = 0;
    } else if(a>b) {
      s = b/a-1;
    }
    s_mean += s;
  }
  return(s_mean/(double)ds_one.n_values);
}

double mean_from_smith_waterman_distance_matrix(unsigned long* matrix,
						dataset ds_one,
						dataset ds_two) {

  size_t i;
  
  double sum = 0.;
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    sum+=(double)matrix[i];
  }
  return(sum/((double)ds_one.n_values*(double)ds_two.n_values));
  
}

double sigma_from_smith_waterman_distance_matrix(unsigned long* matrix,
						 double mean,
						 dataset ds_one,
						 dataset ds_two) {

  size_t i;
  double sigma = 0.;
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    sigma+=((double)matrix[i]-(double)mean)*((double)matrix[i]-(double)mean);
  }
  return(sqrt(sigma)/((double)ds_one.n_values+(double)ds_two.n_values));
}

void shortest_longest_distance_in_matrix(unsigned long* matrix,
					 unsigned long* min,
					 unsigned long* max,
					 dataset ds_one, dataset ds_two) {

  min[0] = ULONG_MAX;
  max[0] = 0;

  size_t i;
  
  for(i=0;i<ds_one.n_values*ds_two.n_values;i++) {
    if(matrix[i] < min[0]) min[0] = matrix[i];
    if(matrix[i] > max[0]) max[0] = matrix[i];
  }
}
				       
