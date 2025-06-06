#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<unistd.h>

#if defined(_SCAN_SMITH_WATERMAN_GPU) || defined(_SCAN_SMITH_WATERMAN_MPI_GPU)
#ifdef __APPLE__
#include<OpenCL/OpenCL.h>
#else
#include<CL/opencl.h>
#endif
#endif

#if defined(_SCAN_SMITH_WATERMAN_MPI_GPU)
#include<mpi.h>
#endif

#if defined(_CLUSTER_KMER_L1) || defined(_CLUSTER_KMER_L2)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#include"dataset.h"
#include"cluster.h"
#include"dbscan.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

void print_arguments() {

  printf("Arguments are: \n"

#if defined (_CLUSTER_KMER_L1) || defined(_CLUSTER_KMER_L2)
	 "   [file] kmer file\n"
#else
	 "   [file] FASTA file containing sequence dataset \n"
#endif
	 "   (float) initial epsilon value \n"
	 "   (float) epsilon increase between cluster search \n"
	 "   (int) minimum numbers of clusters in epsilon neigbourhood \n"
	 "   (string) path/andprefix of split_set_files \n"
	 "   (int) number of threads to use \n"
#if defined (_CLUSTER_PCA)
	 "   (int) dimensions in projection file (PCA obtained from kmers)\n"
	 "   [FILE] PCA projections file obtained i.e. from kers\n"	 
#endif
	 );
}
int main(int argc, char** argv) {

  int i,j,k;
  
  dataset ds;

  float epsilon;
  int minpts;
  
#if defined (_SCAN_SMITH_WATERMAN_GPU) || defined(_SCAN_SMITH_WATERMAN) ||\
  defined(_SCAN_SMITH_WATERMAN_MPI_GPU)
  FILE* fasta_f;
#elif defined (_CLUSTER_PCA)
  FILE* fasta_f;
  FILE* projection_f;
#elif defined (_CLUSTER_KMER_L1) || defined(_CLUSTER_KMER_L2)  
  int kmer_fd;
  FILE* kmer_f;
  data_shape shape;
#endif

#if defined (_SCAN_SMITH_WATERMAN_MPI_GPU)
  int mpi_rank;
  int mpi_size;
  int mpi_kill_switch = -1;
#endif
  
  char split_files_prefix[255];

  float epsilon_start;
  float epsilon_inc;
  int n_threads;
#if defined (_CLUSTER_PCA)
  int dimensions;
#endif

  if(
#if defined (_CLUSTER_PCA)
     argc < 7
#else
     argc < 5
#endif
     ) {
    print_arguments();
    return(1);
  }
  
  sscanf(argv[2], "%f", &epsilon_start);
  sscanf(argv[3], "%f", &epsilon_inc);
  sscanf(argv[4], "%i", &minpts);
  sscanf(argv[5], "%s", split_files_prefix);
  sscanf(argv[6], "%i", &n_threads);
#if defined (_CLUSTER_PCA)
  sscanf(argv[7], "%i", &dimensions);
#endif

#if defined (_SCAN_SMITH_WATERMAN_GPU) || defined(_SCAN_SMITH_WATERMAN)\
  || defined(_SCAN_SMITH_WATERMAN_MPI_GPU)

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

#elif defined (_CLUSTER_PCA)

  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == ( projection_f = fopen(argv[8], "r"))) file_error(argv[8]);
  load_projections_from_file_into_dataset(projection_f,dimensions,&ds);
  fclose(projection_f);

#elif defined (_CLUSTER_KMER_L1) || defined(_CLUSTER_KMER_L2)

  if ( -1 == (kmer_fd = open(argv[1], O_RDONLY))) file_error(argv[1]);
  
  shape = shape_from_kmer_file(kmer_fd);
  kmer_f = fdopen(kmer_fd,"r");
  ds = load_kmer_from_file_into_dataset(kmer_f, shape);
  fclose(kmer_f);

#endif  

  printf("Dataset read!\n");

#if defined(_SCAN_SMITH_WATERMAN_MPI_GPU)
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  if(mpi_rank != 0) {
    adaptive_dbscan_mpi_client(ds, mpi_rank, mpi_size);
    goto finish;
  }
#endif
  
  adaptive_dbscan(
#if defined (_SCAN_SMITH_WATERMAN_GPU)
		   dbscan_SW_GPU,
#elif defined(_SCAN_SMITH_WATERMAN_MPI_GPU)
		   dbscan_SW_GPU_MPI,
#elif defined(_SCAN_SMITH_WATERMAN)
		   dbscan_SW,
#elif defined(_CLUSTER_PCA) || defined(_CLUSTER_KMER_L2)
		   dbscan_L2,
#elif defined(_CLUSTER_KMER_L1)
		   dbscan_L1,
#endif
#if defined(_SCAN_SMITH_WATERMAN_MPI_GPU)
		   mpi_rank,
		   mpi_size,
#endif
		   ds,
		   epsilon_start,
		   epsilon_inc,
		   minpts,
		   split_files_prefix,
		   n_threads);

 finish:
#if defined (_SCAN_SMITH_WATERMAN_MPI_GPU)
  if(mpi_rank == 0) {
    for(i=1;i<mpi_size;i++) {
      MPI_Send(&mpi_kill_switch,1,MPI_INT,i,0,MPI_COMM_WORLD);
    }
  }
#elif defined (_SCAN_SMITH_WATERMAN_GPU) || defined(_SCAN_SMITH_WATERMAN) ||\
  defined (_SCAN_SMITH_WATERMAN_MPI_GPU)
  free_sequences_from_dataset(ds);
#elif defined (_CLUSTER_PCA)
  free_dataset(ds);
#elif defined (_CLUSTER_KMER_L1) || defined(_CLUSTER_KMER_L2)
  free_values_from_dataset(ds);
#endif
#if defined (_SCAN_SMITH_WATERMAN_MPI_GPU)
  MPI_Finalize();
#endif
}
  
  
  

    
    
    
  
    
