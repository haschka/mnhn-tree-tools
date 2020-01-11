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
  cl_device_id* devices;
  cl_platform_id* platform;
  cl_context_properties* contprop;
  cl_program* programs;
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
void opencl_destroy(opencl_stuff ocl);
split_set dbscan_SW_GPU(dataset ds, float epsilon, int minpts,
			opencl_stuff ocl);
void adaptive_dbscan(split_set (*dbscanner) (dataset,
					     float,
					     int,
					     opencl_stuff),
		     dataset ds,
		     float epsilon_start,
		     float epsilon_inc,
		     int minpts,
		     char* split_files_prefix
		     );
#else
void adaptive_dbscan(split_set (*dbscanner) (dataset,
					     float,
					     int),
		     dataset ds,
		     float epsilon_start,
		     float epsilon_inc,
		     int minpts,
		     char* split_files_prefix
		     );
#endif
void create_cluster_files(char* prefix, split_set s, dataset ds);
