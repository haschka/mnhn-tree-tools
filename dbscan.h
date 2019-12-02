/*! \brief Structure storing close or neighbor features of a
 *         feature and used for
 *         the internal workings of dbscan.
 */
typedef struct {
  int * members; /*!< the features in this neighbourhood */
  int n_members; /*!< number of features in this neighborhood */
} neighbors;

split_set dbscan_L1(dataset ds, float epsilon, int minpts);
split_set dbscan_L2(dataset ds, float epsilon, int minpts);
split_set dbscan_SW(dataset ds, float epsilon, int minpts);
void create_cluster_files(char* prefix, split_set s, dataset ds);
