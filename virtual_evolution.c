#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include"dataset.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}


dataset generate_dataset_using_evolution_simulation(int seed,
						    int sequence_length,
						    int n_sequences,
						    int partitions,
						    int max_time,
						    int mutation_rate,
						    int assure_mutation) {
  int i,j;
  dataset ds;
  
  int partition_size;
  
  double random_value_d;
  int random_value_i;

  int partition_count;

  double mutation_time_d;
  int mutation_time_i;

  double mutation_site_d;
  int mutation_site_i;

  double mutation_sequence_d;
  int mutation_sequence_i;

  ds.max_sequence_length = sequence_length;
  ds.n_values = n_sequences;
  
  if( ds.n_values%partitions != 0 ) {
    printf("Size is not eintirely divisible by requested partitions\n");
    _exit(1);
  }
  
  partition_size = ds.n_values/partitions;

  ds.sequences = (char**)malloc(sizeof(char**)*ds.n_values);
  ds.sequence_lengths = (size_t*)malloc(sizeof(size_t)*ds.n_values);
  
  ds.sequences[0] = (char*)malloc(sizeof(char)*(ds.max_sequence_length+1));

  for(i=0;i<ds.n_values;i++) {
    ds.sequence_lengths[i] = (size_t)sequence_length;
  }
  
  srand(seed);

  for(i=0;i<ds.max_sequence_length;i++) {
  
    random_value_d=4*((double)rand()/(double)RAND_MAX);
    random_value_i = (char)random_value_d;

    ds.sequences[0][i] = random_value_i;
    ds.sequences[0][ds.max_sequence_length]= 0;
  }

  for(i=1;i<partition_size;i++) {
    ds.sequences[i] = (char*)malloc(sizeof(char)*(ds.max_sequence_length+1));
    memcpy(ds.sequences[i],ds.sequences[0],sizeof(char)*ds.max_sequence_length);
    ds.sequences[i][ds.max_sequence_length]= 0;
  }

  partition_count = 1;

  while(partition_count < partitions) {

    if (max_time != -1) {
      mutation_time_d = (double)max_time*((double)rand()/(double)RAND_MAX);
    } else {
      mutation_time_d = 1.;
    }
    mutation_time_i = (int)mutation_time_d;

    printf("Partition %3i will be created after %6i cycles of %6i "
	   "mutations/partition\n",
	   partition_count+1, mutation_time_i, mutation_rate);
    
    for(i=0;i<mutation_time_i;i++) {
      for(j=0;j<mutation_rate*partition_count;j++) {
	mutation_site_d =
	  (double)ds.max_sequence_length*(double)(rand()/(double)RAND_MAX);
	mutation_site_i = (int)mutation_site_d;

	mutation_sequence_d =
	  (double)partition_size
	  *(double)partition_count
	  *((double)rand()/(double)RAND_MAX);

	mutation_sequence_i = (int)mutation_sequence_d;
	
	random_value_d=4*((double)rand()/(double)RAND_MAX);
	random_value_i = (char)random_value_d;

	ds.sequences[mutation_sequence_i][mutation_site_i] = random_value_i;

      }
    }

    mutation_sequence_d =
      	  (double)partition_size
	  *(double)partition_count
	  *((double)rand()/(double)RAND_MAX);

    mutation_sequence_i = (int)mutation_sequence_d;

    if(assure_mutation) {
      mutation_site_d =
	(double)ds.max_sequence_length*(double)(rand()/(double)RAND_MAX);
      mutation_site_i = (int)mutation_site_d;
      
      random_value_d=4*((double)rand()/(double)RAND_MAX);
      random_value_i = (char)random_value_d;
      
      ds.sequences[mutation_sequence_i][mutation_site_i] = random_value_i;
    }
    
    for(i=0;i<partition_size;i++) {
      ds.sequences[partition_size*partition_count+i] =
	(char*)malloc(sizeof(char)*(ds.max_sequence_length+1));
      memcpy(ds.sequences[partition_size*partition_count+i],
	     ds.sequences[mutation_sequence_i],
	     sizeof(char)*ds.max_sequence_length);
      ds.sequences[partition_size*partition_count+i][ds.max_sequence_length]=0;
    }
    partition_count++;
  }

  if (max_time != -1) {
    mutation_time_d = (double)max_time*((double)rand()/(double)RAND_MAX);
  } else {
    mutation_time_d = 1.;
  }

  mutation_time_i = (int)mutation_time_d;

  printf("All partitions will be mutated for  %6i cycles of %6i "
	 "mutations/partition\n",mutation_time_i, mutation_rate);
    
  for(i=0;i<mutation_time_i;i++) {
    for(j=0;j<mutation_rate*partition_count;j++) {
      mutation_site_d =
	(double)ds.max_sequence_length*(double)(rand()/(double)RAND_MAX);
      mutation_site_i = (int)mutation_site_d;

      mutation_sequence_d =
	(double)partition_size
	*(double)partition_count
	*((double)rand()/(double)RAND_MAX);

      mutation_sequence_i = (int)mutation_sequence_d;
	
      random_value_d=4*((double)rand()/(double)RAND_MAX);
      random_value_i = (char)random_value_d;

      ds.sequences[mutation_sequence_i][mutation_site_i] = random_value_i;
	
    }
  }

  for(i=0;i<ds.n_values;i++) {
    for(j=0;j<ds.max_sequence_length;j++) {

      switch(ds.sequences[i][j]) {
      case 0:
	ds.sequences[i][j] = 'A';
	break;
      case 1:
	ds.sequences[i][j] = 'C';
	break;
      case 2:
	ds.sequences[i][j] = 'G';
	break;
      case 3:
	ds.sequences[i][j] = 'T';
	break;
      default:
	printf("Fatal error in evolution dataset\n");
	_exit(1);
      }
    }
  }
  return(ds);
}

int main(int argc, char** argv) {

  dataset ds;

  int seed, sequence_length, n_sequences, partitions, max_time, mutation_rate,
    assure_mutation;

  
  FILE* f;
  
  if (argc < 7) {
    printf("Arguments are: \n"
	   " [int] seed\n"
	   " [int] sequence length\n"
	   " [int] number of sequences\n"
	   " [int] number of partitions (amplifications)\n"
	   " [int] maximum time between amplification creation -1 for 1 "
	   " constant.\n"
	   " [int] mutation rate\n"
	   " [int] assure mutation of amplification template\n"
	   " [FASTA] output file\n");
    return(1);
  }

  sscanf(argv[1],"%i",&seed);
  sscanf(argv[2],"%i",&sequence_length);

  sscanf(argv[3],"%i",&n_sequences);
  sscanf(argv[4],"%i",&partitions);
  sscanf(argv[5],"%i",&max_time);
  sscanf(argv[6],"%i",&mutation_rate);
  sscanf(argv[7],"%i",&assure_mutation);
  
  if ( NULL == (f = fopen(argv[8], "w"))) file_error(argv[8]);
  
  ds = generate_dataset_using_evolution_simulation(seed,
						   sequence_length,
						   n_sequences,
						   partitions,
						   max_time,
						   mutation_rate,
						   assure_mutation);

  dataset_to_fasta(f,ds);
  free_sequences_from_dataset(ds);
  fclose(f);

  return(0);
}
