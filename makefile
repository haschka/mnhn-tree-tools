CC=gcc
#CFLAGS=-g -fsanitize=address
#CFLAGS= -g -O1 -march=native -ftree-vectorize
CFLAGS=-g -O2 -march=native -ftree-vectorize
LAPACK=-llapack
MATH=-lm
PTHREAD=-pthread
OPENCL=-lOpenCL
PNG=-lpng

all: fasta2kmer kmer2pca cluster_dbscan_pca cluster_dbscan_kmerL1 \
     cluster_dbscan_kmerL2 cluster_dbscan_SW cluster_dbscan_SW_GPU \
     compareSW silhouette consens sequence_multiplicity adaptive_clustering_SW \
     adaptive_clustering_SW_GPU adaptive_clustering_PCA \
     adaptive_clustering_kmer_L1 adaptive_clustering_kmer_L2 \
     split_set_to_fasta print_connections \
     split_set_to_matrix_line split_set_to_matrix_annotation \
     split_sets_to_newick virtual_evolution simulation_verification \
     find_sequence_in_split_sets tree_map_for_sequence pca2densitymap

compare.o: compare.c compare.h dataset.h smith-waterman.h
	$(CC) $(CFLAGS) -c comparison.c -o comparison
smith_waterman.o: smith-waterman.c smith-waterman.h
	$(CC) $(CFLAGS) -c smith-waterman.c -o smith_waterman.o
binary_array.o: binary_array.c binary_array.h
	$(CC) $(CFLAGS) -c binary_array.c -o binary_array.o
cluster_io.o: cluster_io.c cluster.h dataset.h binary_array.h
	$(CC) $(CFLAGS) -c cluster_io.c -o cluster_io.o
dbscan_SW.o: dbscan.c dbscan.h cluster.h binary_array.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_SMITH_WATERMAN -o dbscan_SW.o
dbscan_SW_GPU.o: dbscan.c dbscan.h cluster.h binary_array.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_SMITH_WATERMAN_GPU \
 -o dbscan_SW_GPU.o
dbscan_L1.o: dbscan.c dbscan.h cluster.h binary_array.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_L1 -o dbscan_L1.o
dbscan_L2.o: dbscan.c dbscan.h cluster.h binary_array.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_L2 -o dbscan_L2.o 
dataset.o: dataset.c dataset.h binary_array.h
	$(CC) $(CFLAGS) -c dataset.c -o dataset.o
kmers.o: kmers.c dataset.h kmers.h
	$(CC) $(CFLAGS) -c kmers.c -o kmers.o
density.o: density.c density.h dataset.h
	$(CC) $(CFLAGS) -c density.c -o density.o

fasta2kmer: fasta2kmer.c dataset.h kmers.h dataset.o kmers.o binary_array.o
	$(CC) $(CFLAGS) fasta2kmer.c -o ./bin/fasta2kmer dataset.o kmers.o \
 binary_array.o $(PTHREAD) $(MATH)

kmer2pca: kmer2pca.c 
	$(CC) $(CFLAGS) kmer2pca.c -o ./bin/kmer2pca -mavx $(PTHREAD) $(MATH) \
 $(LAPACK)

pca2densitymap: pca2densitymap.c dataset.h binary_array.h density.h dataset.o \
                binary_array.o density.o 
	$(CC) $(CFLAGS) pca2densitymap.c -o ./bin/pca2densitymap \
 dataset.o binary_array.o density.o $(PNG) $(MATH)

virtual_evolution: virtual_evolution.c dataset.h binary_array.h \
                   dataset.o binary_array.o
	$(CC) $(CFLAGS) virtual_evolution.c -o ./bin/virtual_evolution \
 dataset.o binary_array.o $(MATH)

simulation_verification: simulation_verification.c dataset.h binary_array.h \
                         cluster.h dataset.o binary_array.o cluster_io.o
	$(CC) $(CFLAGS) simulation_verification.c \
 -o ./bin/simulation_verification dataset.o binary_array.o cluster_io.o $(MATH)

find_sequence_in_split_sets: find_sequence_in_split_sets.c dataset.h \
                             binary_array.h cluster.h dataset.o binary_array.o \
                             cluster_io.o
	$(CC) $(CFLAGS) find_sequence_in_split_sets.c \
 -o ./bin/find_sequence_in_split_sets dataset.o binary_array.o cluster_io.o \
 $(MATH)

tree_map_for_sequence: tree_map_for_sequence.c dataset.h \
                       binary_array.h cluster.h dataset.o binary_array.o \
                       cluster_io.o
	$(CC) $(CFLAGS) tree_map_for_sequence.c \
 -o ./bin/tree_map_for_sequence dataset.o binary_array.o cluster_io.o \
 $(MATH)


cluster_dbscan_pca: cluster_dbscan_pca.c dbscan.h dataset.h cluster.h \
                    dataset.o cluster_io.o dbscan_L2.o binary_array.o
	$(CC) $(CFLAGS) cluster_dbscan_pca.c -o ./bin/cluster_dbscan_pca \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o $(MATH)

cluster_dbscan_kmerL1: cluster_dbscan_kmer.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_L1.o binary_array.o
	$(CC) $(CFLAGS) cluster_dbscan_kmer.c -o ./bin/cluster_dbscan_kmerL1 \
 dataset.o cluster_io.o dbscan_L1.o binary_array.o -D_CLUSTER_KMER_L1 $(MATH)

cluster_dbscan_kmerL2: cluster_dbscan_kmer.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_L2.o binary_array.o
	$(CC) $(CFLAGS) cluster_dbscan_kmer.c -o ./bin/cluster_dbscan_kmerL2 \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o -D_CLUSTER_KMER_L2 $(MATH)

cluster_dbscan_SW: cluster_dbscan_SW.c dbscan.h dataset.h cluster.h \
                   dataset.o cluster_io.o dbscan_SW.o binary_array.o \
                   smith_waterman.o
	$(CC) $(CFLAGS) cluster_dbscan_SW.c -o ./bin/cluster_dbscan_SW \
 dataset.o cluster_io.o dbscan_SW.o binary_array.o smith_waterman.o $(MATH)

cluster_dbscan_SW_GPU: cluster_dbscan_SW.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_SW_GPU.o binary_array.o
	$(CC) $(CFLAGS) cluster_dbscan_SW.c -o ./bin/cluster_dbscan_SW_GPU \
 dataset.o cluster_io.o dbscan_SW_GPU.o binary_array.o $(MATH) $(OPENCL) \
 -D_SCAN_SMITH_WATERMAN_GPU

adaptive_clustering_SW: adaptive_clustering.c dbscan.h dataset.h cluster.h \
                     dataset.o cluster_io.o dbscan_SW.o binary_array.o \
                     smith_waterman.o
	$(CC) $(CFLAGS) adaptive_clustering.c -o ./bin/adaptive_clustering_SW \
 dataset.o cluster_io.o dbscan_SW.o binary_array.o smith_waterman.o $(MATH) \
 -D_SCAN_SMITH_WATERMAN

adaptive_clustering_SW_GPU: adaptive_clustering.c dbscan.h dataset.h cluster.h \
	                 dataset.o cluster_io.o dbscan_SW_GPU.o binary_array.o \
                         smith_waterman.o
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_SW_GPU \
 dataset.o cluster_io.o dbscan_SW_GPU.o binary_array.o smith_waterman.o \
 $(MATH) $(OPENCL) -D_SCAN_SMITH_WATERMAN_GPU

adaptive_clustering_PCA: adaptive_clustering.c dbscan.h dataset.h cluster.h \
	                 dataset.o cluster_io.o dbscan_L2.o binary_array.o \
                         smith_waterman.o
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_PCA \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o smith_waterman.o \
 $(MATH) -D_CLUSTER_PCA

adaptive_clustering_kmer_L2: adaptive_clustering.c dbscan.h dataset.h cluster.h\
	                 dataset.o cluster_io.o dbscan_L2.o binary_array.o \
                         smith_waterman.o
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_kmer_L2 \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o smith_waterman.o \
 $(MATH) -D_CLUSTER_KMER_L2

adaptive_clustering_kmer_L1: adaptive_clustering.c dbscan.h dataset.h cluster.h\
	                 dataset.o cluster_io.o dbscan_L1.o binary_array.o \
                         smith_waterman.o
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_kmer_L1 \
 dataset.o cluster_io.o dbscan_L1.o binary_array.o smith_waterman.o \
 $(MATH) -D_CLUSTER_KMER_L1

split_set_to_fasta: split_set_to_fasta.c dataset.h cluster.h dataset.o \
                    cluster_io.o binary_array.o
	$(CC) $(CFLAGS) split_set_to_fasta.c -o ./bin/split_set_to_fasta \
 dataset.o cluster_io.o binary_array.o $(MATH)

split_sets_to_newick: split_sets_to_newick.c dataset.h cluster.h dataset.o \
                    cluster_io.o binary_array.o
	$(CC) $(CFLAGS) split_sets_to_newick.c -o ./bin/split_sets_to_newick \
 dataset.o cluster_io.o binary_array.o $(MATH)

split_set_to_matrix_line: split_set_to_matrix_line.c dataset.h cluster.h \
                          dataset.o cluster_io.o binary_array.o
	$(CC) $(CFLAGS) split_set_to_matrix_line.c \
 -o ./bin/split_set_to_matrix_line dataset.o cluster_io.o binary_array.o $(MATH)

split_set_to_matrix_annotation: split_set_to_matrix_line.c dataset.h cluster.h \
                                binary_array.o cluster_io.o dataset.o
	$(CC) $(CFLAGS) split_set_to_matrix_annotation.c \
 -o ./bin/split_set_to_matrix_annotation dataset.o cluster_io.o \
 binary_array.o $(MATH)

print_connections: print_connections.c dataset.h cluster.h dataset.o \
                    cluster_io.o binary_array.o
	$(CC) $(CFLAGS) print_connections.c -o ./bin/print_connections \
 dataset.o cluster_io.o binary_array.o $(MATH)

compareSW: compareSW.c dataset.h comparison.h smith_waterman.o comparison.o \
           dataset.o binary_array.o
	$(CC) $(CFLAGS) compareSW.c -o ./bin/compareSW smith_waterman.o \
 comparison.o dataset.o binary_array.o $(MATH) $(PTHREAD)

silhouette: dataset.h comparison.h dataset.h dataset.o comparison.o \
            smith_waterman.o binary_array.o
	$(CC) $(CFLAGS) silhouette.c -o ./bin/silhouette smith_waterman.o \
 comparison.o dataset.o binary_array.o $(MATH) $(PTHREAD)

consens: consens.c dataset.h dataset.o binary_array.o
	$(CC) $(CFLAGS) consens.c -o ./bin/consens dataset.o binary_array.o \
 $(MATH)

sequence_multiplicity: sequence_multiplicity.c dataset.h binary_array.o \
                       dataset.o
	$(CC) $(CFLAGS) sequence_multiplicity.c -o ./bin/sequence_multiplicity \
 binary_array.o dataset.o $(MATH)

clean:
	rm ./bin/* *.o
