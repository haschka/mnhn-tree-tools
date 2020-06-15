#include<math.h>

double vol_hypersphere(int dimension, double epsilon) {

  double dim_f = (double) dimension;
  
  return( pow(epsilon,dim_f)*((pow(M_PI,dim_f*0.5)) / (gamma(dim_f*0.5+1))) );
}

double vol_hypercube(int dimension, double epsilon) {

  int i;
  double vol = 1.;
  for(i=0;i<dimension;i++) {
    vol *= epsilon;
  }
  return(vol);
  
}
  
