/*! \brief Structure for a n dimenstional dataset composed of
 *         non missing floating point values.
 */
typedef struct {
  char** sequences;
  size_t* sequence_lengths;
  void** binary_sequences;
  int n_values; /*!< number of samples in the dataset*/
  float** values; /*!< supplimental data to sequencies, i.e. used for kmers */
  int n_dimensions; /*!< dimensions of supplimental data*/
  size_t max_sequence_length;
} dataset;

typedef struct {
  size_t n_features;
  size_t n_samples;
} data_shape;

dataset dataset_from_fasta(FILE* in);
void free_dataset(dataset ds);
void free_values_from_dataset(dataset ds);
void free_sequences_from_dataset(dataset ds);
void load_projections_from_file_into_dataset(FILE* projections,
					     size_t dimensions,
					     dataset* ds);

dataset load_kmer_from_file_into_dataset(FILE* in_file, data_shape shape);
data_shape shape_from_kmer_file(int infile);
