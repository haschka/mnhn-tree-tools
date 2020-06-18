/*! \file cluster.h
 *  \brief defines functions for cluster I/O and basic cluster manipulation and
 *         relation. 
 */ 

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

/*! \brief A structure used in the generation of trees
 *  \var id A string holding the name of this threenode
 *  \var length the length of the branch to this node (parent to this)
 *  \var n_members The number of sequencies ( objects ) this node holds
 *  \var child a pointer to a node on a lower level this one.
 *  \var parent a pointer to the parent node of this one
 *  \var neighbor a pointer to a neiboring node to this one.
 */
typedef struct tree_node{
  char* id;
  double length;
  int n_members;
  struct tree_node* child;
  struct tree_node* neighbor;
  struct tree_node* parent;
} tree_node;

/*! \brief creates files in fasta format for each cluster indexed in 
 *         a given split_set.
 *  This function allows to print out fasta files from clusters obtained 
 *  in a clustering run. For this function to work the dataset given 
 *  has to hold the sequences of the clusters. Multiple files are generated: 
 *  one per cluster. Each file contains the name prefixN where N is the number
 *  of the cluster as it occurs in the split_set. N is prined out in the printf
 *  format "%04d" and hence contains leading zeros. 
 *  \param prefix a file path and prefix for the files to be generated. 
 *                for each cluster a file with the name prefixN will be
 *                generated.
 *                where N is the number of the cluster as it occurs in the 
 *                split_set s
 *  \param s the split_set containing the clusters to printed out as 
 *           individual fasta files.
 *  \param ds the dataset that the clusters / split_set was obtained from
 */ 
void create_cluster_files(char* prefix, split_set s, dataset ds);

/*! \brief Creates cluster files containing not the sequence representation
 *         but a different value representation for each sequence, or datapoint.
 *  This function basically workes just as the creates cluster files function, 
 *  but does not generate files in fasta format but prints out the values 
 *  associated to each sequence / datapoint. For instance a pca projection of 
 *  a kmer representation of the sequence in question. Each data point is
 *  on a new line in the data output file. Each values for this data point 
 *  ( sequence ) is tabulator seperated. For this function to work the 
 *  dataset provided to this function must hold the values to be printed out. 
 *  \param prefix a file path and prefix for the files to be generated. 
 *                for each cluster a file with the name prefixN will be
 *                generated.
 *                where N is the number of the cluster as it occurs in the 
 *                split_set s
 *  \param s the split_set containing the clusters to printed out as 
 *           individual values files.
 *  \param ds the dataset that the clusters / split_set was obtained from. 
 *            The dataset has to hold the values to be printed out.
 */ 
void create_cluster_files_with_values(char* prefix, split_set s, dataset ds);

void create_single_cluster_file(char* filename, cluster cl, dataset ds);
void create_single_cluster_file_with_values(char* filename, cluster cl,
					    dataset ds);

/*! \brief Creates a cluster from the intersection of two clusters.
 *  \param a A cluster.
 *  \param b An other cluster.
 *  \return the intersection ( the datapoints / sequences contained in 
 *          cluster a and cluster b.
 */
cluster intersection_of_clusters(cluster a, cluster b);

/*! \brief Creates a cluster containing all the data not contained in the 
 *         split_set given. 
 *  This function provides a cluster of the data not covered by a set of 
 *  clusters. I.e. the outliers of a dbscan run. 
 *  \param s A split_set.
 *  \param ds A dataset where the split_set s operates on
 *  \return A cluster containing all the values / sequences not 
 *          covered in the split_set s. 
 */
cluster data_not_in_clusters(split_set s, dataset ds);

/*! \brief Permits to obtain related clusters between to split sets.
 *  This function finds out if a cluster contained in one split set 
 *  is contained in the other split set and if yes provides an index 
 *  to trace this relation. In the current implementation a cluster from 
 *  one split set to the other split set is related if it contains 80% or 
 *  more of the datapoints / sequences of a cluster in a different split set. 
 *  ancient clusters are searched in new clusters. From a runtime perspective, 
 *  the splitset containing less clusters should be the new one. 
 *  cluser_connections* is organized that for each cluster in the new split set
 *  a cluster_connection is generated holding the connections from the ancient
 *  to the new split set. 
 *  \param ancient a split set
 *  \param new an other split set preferably containing less clusters than
 *             ancient
 *  \return A pointer to an array of cluster_connections, holding one 
 *          cluster_connections structure for each cluster in the new 
 *          split set. 
 */
cluster_connections* generate_split_set_relation(split_set ancient,
						 split_set new);


/*! \brief permits to create a binary file from a split set and 
 *         all the cluster files it containes on non volatile disk storage. 
 *  \param filename The path/filename where to store the split_set to. 
 *  \param s The split set to be stored. 
 */
void store_split_set(char* filename, split_set s);

/*! \brief recover a stored split set and its clusters from a disk.
 *  \param the file name of the stored split_set. 
 *  \return the recovered split_set.
 */
split_set read_split_set(char* filename);

split_set filtered_split_set_by_min_size(split_set s_in, int min_size);

void print_cluster_matrix_view_annotation(FILE* f, dataset ds);
void print_cluster_matrix_view(FILE*f, split_set s, dataset ds);

/* \brief a function that frees a split set and all clusters contained indexed
 *        by the split set. 
 * \param s the split_set to be freed.
 */
void free_split_set_and_associated_clusters(split_set s);

/* \brief a function to generate trees from split_set s for instance
 *        obtained from an adaptive clustering run.
 * This basically provides a pointer linked recursive tree, which is useful
 * for applications like creating a tree in newick format, from previously
 * caluculated connections (cluster correspondances) between splits_sets.
 * The root node (sets[n-1]) if n sets are provided. Each set below the root
 * node contains a pointer to its parent node. Even if a node has k children
 * only the childnode[k-1] is indicated as child. If multiple childrens exist, 
 * hence the other childrens are accessible using the neigbour pointer.
 * For each child_node[i] besides child_node[0] its neighbor is 
 * child_node[i-1]. child_node[0] contains no neighbor and the neighbor pointer
 * is set to NULL. 
 * \param n_layers the number of split_sets, layers that this tree is 
 *                 built from
 * \param c the connections between all the layers split_sets.
 *          Such connections can be obtained by calling
 *          generate_split_set_relation. If there are n split_sets
 *          then c has to hold n-1 (cluster_connections*), where
 *          connection[i] = generate_split_set_relation(sets[i],sets[i+1]);
 * \param sets the split_sets obtained from a clusting run or adaptive 
 *             clustering run.
 * \param tree_lengths an array to set the branch lengths of the trees 
 *                     between the layers. If tree_lengths == NULL 
 *                     The branch lengths are 1 everywhere.
 * \return A tree node structure holding the entire through recursion
 */
tree_node* generate_tree(int n_layers, cluster_connections** c,
			 split_set* sets, double* tree_lengths);

/* \brief prints a tree in newick format to a given file
 * This function uses a recursive tree traversal method which is optimal
 * in order to generate newick files
 * \param f an opened writeable file to write to
 * \param root the tree to write out, i.e. obtained by generate_tree
 */
void print_tree(FILE* f, tree_node* root);

double* get_epsilon_dist_from_adaptive_clustering_output(FILE* f,
							 int n_layers);

double* densities_from_epsilons(int type, int dimensions,
				double* epsilons, int n_epsilons,
				int min_points);


/* \brief creates a cluster that holds, provides an index for identical 
 *        sequences.
 * \param ds the dataset to create this cluster from.
 * \param seq the sequence that is searched in the dataset, and if it occures
 *            its index is added to cluster. If it occurs multiple times,
 *            all these sequences are added to the cluster
 * \param seq_len the length of the sequence pointed to by seq. 
 * \return a cluster structure containing all sequences that are identical to
 *         the sequence in seq.
 */
cluster cluster_from_sequence_in_dataset(dataset ds,
					 char* seq,
					 size_t seq_len);


double* array_deltas(double* array, int array_length);
