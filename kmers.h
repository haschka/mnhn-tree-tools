/*! \file kmers.h
 *  \brief contains functions to work with kmers and a structure 
 *         holding kmer frequencies for a dataset containing sequences
 */


/*! \brief a structure holding the frequencies of different kmers in 
 *           sequencies
 *  \var frequencies the frequencies table for each sequence and kmer
 *                   organized in frequencies[0 .. n_seq-1][0 .. n_kmers-1],
 *                   where the first index is the sequence index and the
 *                   second is the k_mer index
 *  \var n_kmers the number of k_mers used to obtain this frequency table
 *  \var n_seq the number of sequences used to obtain this frequency table
 */
typedef struct {
  short int** frequencies;
  size_t n_kmers;
  size_t n_seq;
} kmer_frequencies;

/*! \brief a function that returns the number of possible k_mers for 
 *         a specific k, hence k-mer length.
 *  \param kmer_length the k length paramerter.
 *  \return the number of possible k_mers. 
 */
size_t number_of_kmers_from_length(size_t kmer_length);

/*! \brief a function that generates all possible kmers for a given length.
 *         Warning this function only works for k_mer lenghts up to 16.
 *  \param kmer_length the length of the k_mers to be obtained
 *  \return DNA strings of the k_mers
 */
char** gen_dna_kmers(size_t kmer_length);

/*! \brief a function that returns the frequency for a given kmer in a 
 *         given sequence.
 *  \param kmer the DNA string of the k_mer in question
 *  \param sequence the DNA string of the Sequence the k_mer should be found on.
 *  \param kmer_length the length of the k_mer in question
 *  \return the frequency that the k_mer was found in the sequence.
 */
short int frequence_of_kmer_in_sequence(char* kmer, char* sequence,
					size_t kmer_length);

/*! \brief a function to free the memory that has been allocated for k_mers.
 *  \param kmers a pointer to the k_mers in question
 *  \param kmer_length the length of the k_mers in question
 */
void free_kmers(char** kmers, size_t kmer_length);

void* frequencies_from_dataset_thread_handler(void* arg);

/*! \brief Allows you to obtain the frequencies of k_mers occuring in 
 *         sequences accross a dataset.
 *  This function generates a kmer_frequencies structure that corresponds
 *  to a k_mer representation of the given dataset.
 *  \param ds the dataset to containing sequences to represent as 
 *            k_mer frequency vectors
 *  \param kmer_length the length of the k_mers to obtain the frequencies 
 *                     form. Currently only kmers up to the length of 15 are
 *                     supported.
 *  \param n_threads the number of threads to use for this run
 */
kmer_frequencies frequencies_from_dataset(dataset ds, size_t kmer_length,
					  int n_threads);

/*! \brief Write kmer_frequencies out to a file. 
 *  This function allows you to write the frequencies obtained for instance
 *  from a dataset of sequences out to a file
 *  \param f A file pointer that is opened and writeable.
 *  \param freq A kmer_frequencies structure holding the frequencies 
 *              to be written to the file.
 */
void write_kmer_base(FILE* f,kmer_frequencies freq);

/*! \brief a function that frees the memory allocated to a kmer_frequencies
 *         structure.
 *  \param freq the structure to be freed
 */
void free_kmer_frequencies(kmer_frequencies freq);

