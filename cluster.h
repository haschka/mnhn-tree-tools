/*! \brief Structure storing a single cluster to be found by
 *         clustering algorithms
 *         such as dbscan.
 */
typedef struct {
  int * members;  /*!< integer containing the number of each sample that is
		   *   part of this cluster */
  int n_members;  /*!< number of samples represented in this cluster */
  int id;         /*!< unique cluster id in a set of clusters */
} cluster;

/*! \brief Structure storing multiple clusters of a dataset to be found by
 *         clustering algorithms such as dbscan.
 */
typedef struct {
  int n_clusters; /*!< number of clusters to be contained in this set */
  cluster* clusters; /*!< the clusters contained within this set */
} split_set;

/*! \brief Structure for storing connections between clusters, i.e. to build
 *         dodendograms
 */
typedef struct {
  int* connections;  /*!< array holding the connections for a single cluster */
  int n_connections; /*!< the number of connections for a single cluster */
} cluster_connections;

typedef struct tree_node{
  char* id;
  double length;
  int n_members;
  struct tree_node* child;
  struct tree_node* neighbor;
  struct tree_node* parent;
} tree_node;


void create_cluster_files(char* prefix, split_set s, dataset ds);
void create_cluster_files_with_values(char* prefix, split_set s, dataset ds);
void create_single_cluster_file(char* filename, cluster cl, dataset ds);
void create_single_cluster_file_with_values(char* filename, cluster cl,
					    dataset ds);
cluster intersection_of_clusters(cluster a, cluster b);
cluster data_not_in_clusters(split_set s, dataset ds);
cluster_connections* generate_split_set_relation(split_set ancient,
						 split_set new);


void store_split_set(char* filename, split_set s);
split_set read_split_set(char* filename);
void print_cluster_matrix_view_annotation(FILE* f, dataset ds);
void print_cluster_matrix_view(FILE*f, split_set s, dataset ds);

void free_split_set_and_associated_clusters(split_set s);
tree_node* generate_tree(int n_layers, cluster_connections** c,
			 split_set* sets);
void print_tree(FILE* f, tree_node* root);
cluster cluster_from_sequence_in_dataset(dataset ds,
					 char* seq,
					 size_t seq_len);
