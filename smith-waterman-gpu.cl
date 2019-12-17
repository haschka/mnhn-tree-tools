int SCORE=4;
int GAP=3; 
    
int substitute(char a, char b) {

  if(a == b) {
    return(SCORE);
  } else {
    return (-SCORE);
  }
}

__kernel void gpuwaterman(__global char* sequences,
			  __global size_t* lengths,
			  __global int* ret_vals,
			  int sequence_offset,
			  int idx_current_seq ) {


  __private int matrix[600];
  __private int comp[3];
  
  int i,j,k;

  int free_matrix = 0;

  int idx = get_global_id(0);
  
  __global char* A = sequences+idx*sequence_offset;
  __global char* B = sequences+idx_current_seq*sequence_offset;

  //  printf("%i \n",);
  printf("%s \n",A);
  printf("%s \n",B);
  printf("i %i \n", idx);
  
  int len_a = lengths[idx];
  int len_b = lengths[idx_current_seq];
  
  int rank_a = len_a+1;
  int rank_b = len_b+1;

  int optimum, max_max = 0;
  int max;

  // printf("a %i\n", len_a);
  //printf("b %i\n", len_b);
  
  if ( len_a > len_b ) {
    optimum = len_a*SCORE;
  } else {
    optimum = len_b*SCORE;
  }

  for(i=0;i<rank_a;i++) {
    matrix[i] = 0;
  }

  for(j=0;j<3;j++) {
    matrix[j*rank_a] = 0;
  }

  for(j=1;j<rank_b;j++) {
    for(i=1;i<rank_a;i++) {
      comp[0] = matrix[(i-1)]+substitute(A[i],B[j]);
      comp[1] = matrix[rank_a+(i-1)]-GAP;
      comp[2] = matrix[i]-GAP;
      //printf("1: %i\n",matrix[i-1]);
      //printf("2: %i\n",matrix[rank_a]);
      //printf("3: %i\n",comp[2]);

      
      max = 0;

      for(k=0;k<3;k++) {
	if(max < comp[k]) {
	  max = comp[k];
	}
      }

      matrix[rank_a+i] = max;
      //  printf("max %i\n",max); 
      if(max_max < max) {
	max_max = max;
      }
    }
    for(i=0;i<rank_a;i++) {
      matrix[i] = matrix[rank_a+i];
    }
  }
  //  printf("m %i\n", max_max);
  //printf("o %i\n", optimum);
    
  ret_vals[idx] = optimum - max_max;
}
      
		       

  
  
  

			 
		       
