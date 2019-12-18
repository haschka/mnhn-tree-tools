/*! \brief Structure storing close or neighbor features of a
 *         feature and used for
 *         the internal workings of dbscan.
 */
typedef struct {
  int * members; /*!< the features in this neighbourhood */
  int n_members; /*!< number of features in this neighborhood */
} neighbors;

#if defined(_SCAN_SMITH_WATERMAN_GPU)
typedef struct {
  cl_command_queue* cmdq;
  cl_kernel* kernel;
  cl_context* contexts;
  cl_mem* acc_sequences;
  cl_mem* acc_sequence_lengths;
  cl_mem* acc_distances;
  int * local_distances;
  int num_devs;
} opencl_stuff;
#endif


split_set dbscan_L1(dataset ds, float epsilon, int minpts);
split_set dbscan_L2(dataset ds, float epsilon, int minpts);
split_set dbscan_SW(dataset ds, float epsilon, int minpts);
#if defined(_SCAN_SMITH_WATERMAN_GPU)
opencl_stuff opencl_initialization(dataset ds);
split_set dbscan_SW_GPU(dataset ds, float epsilon, int minpts,
			opencl_stuff ocl);
#endif
void create_cluster_files(char* prefix, split_set s, dataset ds);
