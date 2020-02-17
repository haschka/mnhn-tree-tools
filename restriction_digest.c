#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"dataset.h"

dataset digest_XbaI(char * sequence, size_t sequence_length) {

  size_t i;
  size_t* sites = (size_t*)malloc(sizeof(size_t)*1000);

  size_t n_sites = 0;

  size_t max_size;
  
  dataset ds;
  
  for(i=1;i<sequence_length;i++) {

    if(sequence[i-1] == 'T' &&
       sequence[  i] == 'C' &&
       sequence[i+1] == 'T' &&
       sequence[i+2] == 'A' &&
       sequence[i+3] == 'G' &&
       sequence[i+4] == 'A') {

      n_sites++;
      if(n_sites %1000 == 0) {
	sites = (size_t*)realloc(sites,sizeof(size_t)*(n_sites+1000));
      }
      sites[n_sites-1] = i;
    }
  }


  ds.n_values = n_sites+1;
  ds.sequences = (char**)malloc(sizeof(char*)*(ds.n_values));
  ds.sequence_lengths = (size_t*)malloc(sizeof(size_t)*(ds.n_values));
    
  
  if(n_sites == 0) {
    ds.sequences[0] = (char*)malloc(sizeof(char)*(sequence_length+1));
    memcpy(ds.sequences[0],sequence,sequence_length);
    ds.sequences[0][sequence_length] = 0;
    ds.sequence_lengths[0] = sequence_length;
  }

  if(n_sites > 0) {
    
    /*first*/
    ds.sequences[0] = (char*)malloc(sizeof(char)*(sites[0])+1);
    memcpy(ds.sequences[0],sequence,sizeof(char)*(sites[0]));
    ds.sequences[0][sites[0]] = 0;
    ds.sequence_lengths[0] = sites[0];
    /*last*/
    ds.sequence_lengths[n_sites] = sequence_length-sites[n_sites-1];
    ds.sequences[n_sites] =
      (char*)malloc(sizeof(char)*(ds.sequence_lengths[n_sites]+1));
    memcpy(ds.sequences[n_sites],sequence+sites[n_sites-1],sizeof(char)*
	   ds.sequence_lengths[n_sites]);		    
    /*inbetween*/
    for(i=0;i<n_sites-1;i++) {
      ds.sequence_lengths[i+1] = sites[i+1]-sites[i];
      ds.sequences[i+1] = (char*)malloc(sizeof(char)*
					ds.sequence_lengths[i+1]+1);
      
      memcpy(ds.sequences[i+1],sequence+sites[i],ds.sequence_lengths[i+1]);
      ds.sequences[i+1][ds.sequence_lengths[i+1]] = 0;
    }
  }

  max_size = 0;
  for(i=0;i<ds.n_values;i++) {
    if(ds.sequence_lengths[i] > max_size) max_size = ds.sequence_lengths[i];
  }
  ds.max_sequence_length = max_size;

  return(ds);
}
	
    
