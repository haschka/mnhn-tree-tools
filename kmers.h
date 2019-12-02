typedef struct {
  short int** frequencies;
  size_t n_kmers;
  size_t n_seq;
} kmer_frequencies;

size_t number_of_kmers_from_length(size_t kmer_length);

char** gen_dna_kmers(size_t kmer_length);

short int frequence_of_kmer_in_sequence(char* kmer, char* sequence,
					size_t kmer_length);

void free_kmers(char** kmers, size_t kmer_length);

void* frequencies_from_dataset_thread_handler(void* arg);

kmer_frequencies frequencies_from_dataset(dataset ds, size_t kmer_length,
					  int n_threads);

void write_kmer_base(FILE* f,kmer_frequencies freq);

void free_kmer_frequencies(kmer_frequencies freq);

