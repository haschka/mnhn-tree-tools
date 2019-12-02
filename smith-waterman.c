#include <string.h>
#include <stdlib.h>

#define SCORE 4
#define GAP 3

static inline int substitute(char a, char b) {

  if(a == b) {
    return(SCORE);
  } else {
    return (-SCORE);
  }
}

unsigned int score(char* A, char* B, size_t len_a, size_t len_b) {

  int i,j,k;
  
  size_t rank_a = len_a+1;
  size_t rank_b = len_b+1;
  
  int* matrix = (int*)malloc(sizeof(int)*(rank_a)*(rank_b));

  int comp[3];

  int max;

  int max_max = 0;

  int optimum;

  if (len_a > len_b) {
    optimum = len_a*SCORE;
  } else {
    optimum = len_b*SCORE;
  }
  
  for(i=0;i<rank_a;i++) {
    matrix[i] = 0;
  }

  for(j=0;j<rank_b;j++) {
    matrix[j*rank_a] = 0;
  }

  for(j=1;j<rank_b;j++) {
    for(i=1;i<rank_a;i++) {
      comp[0] = matrix[(j-1)*rank_a+(i-1)]+substitute(A[i],B[j]);
      comp[1] = matrix[j*rank_a+(i-1)]-GAP;
      comp[2] = matrix[(j-1)*rank_a+i]-GAP;

      max = 0;

      for(k=0;k<3;k++) {
	if (max < comp[k]) {
	  max = comp[k];
	}
      }

      matrix[j*rank_a+i] = max;
      
      if(max_max < max) {
	max_max = max;
      }
    }
  }

  free(matrix);
  return(optimum-max_max);
}



  
  
