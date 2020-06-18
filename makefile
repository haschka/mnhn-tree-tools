CC=gcc
MPICC=mpicc
CFLAGS=-g
#CFLAGS=-g -fsanitize=address 
#CFLAGS= -g -O1 -march=native -ftree-vectorize
#CFLAGS= -O3 -mcpu=750 -mtune=750 -fomit-frame-pointer 
#CFLAGS= -O2 -march=native -ftree-vectorize -ftree-loop-linear -fomit-frame-pointer
LAPACK=-llapack
MATH=-lm
PTHREAD=-pthread
OPENCL=-lOpenCL
PNG=-lpng
MPI=-lmpi
SDL=$(shell sdl2-config --cflags) $(shell sdl2-config --libs)

all: fasta2kmer kmer2pca cluster_dbscan_pca cluster_dbscan_kmerL1 \
     cluster_dbscan_kmerL2 cluster_dbscan_SW cluster_dbscan_SW_GPU \
     compareSW silhouette compare_norms consens sequence_multiplicity \
     adaptive_clustering_SW \
     adaptive_clustering_SW_GPU adaptive_clustering_PCA \
     adaptive_clustering_kmer_L1 adaptive_clustering_kmer_L2 \
     adaptive_clustering_SW_MPI_GPU \
     split_set_to_fasta print_connections \
     split_set_to_matrix_line split_set_to_matrix_annotation \
     split_sets_to_newick split_set_to_projections \
     virtual_evolution virtual_evolution_controlled simulation_verification \
     find_sequence_in_split_sets tree_map_for_sequence tree_map_for_split_set \
     filter_split_sets_by_min \
     pca2densitymap pca2densityfile \
     reverse_with_mask reverse_complement_with_mask \
     find_closest_sequence_SW split_set_from_annotation \
     pca_visual_extract digest_XbaI digest_XmnI digest_HindIII \
     find_satellite

comparison.o: comparison.c comparison.h dataset.h smith-waterman.h
	$(CC) $(CFLAGS) -c comparison.c -o comparison.o
smith_waterman.o: smith-waterman.c smith-waterman.h
	$(CC) $(CFLAGS) -c smith-waterman.c -o smith_waterman.o
binary_array.o: binary_array.c binary_array.h
	$(CC) $(CFLAGS) -c binary_array.c -o binary_array.o
volumes.o: volumes.c volumes.h
	$(CC) $(CFLAGS) -c volumes.c -o volumes.o
cluster_io.o: cluster_io.c cluster.h dataset.h binary_array.h volumes.h
	$(CC) $(CFLAGS) -c cluster_io.c -o cluster_io.o
dbscan_SW.o: dbscan.c dbscan.h cluster.h binary_array.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_SMITH_WATERMAN -o dbscan_SW.o
dbscan_SW_GPU.o: dbscan.c dbscan.h cluster.h binary_array.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_SMITH_WATERMAN_GPU \
 -o dbscan_SW_GPU.o
dbscan_SW_MPI_GPU.o: dbscan.c dbscan.h cluster.h binary_array.h
	$(MPICC) $(CFLAGS) -c dbscan.c -D_SCAN_SMITH_WATERMAN_MPI_GPU \
 -o dbscan_SW_MPI_GPU.o
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
filter.o: filter.c filter.h dataset.h smith-waterman.h
	$(CC) $(CFLAGS) -c filter.c -o filter.o
restriction_digest_XbaI.o: restriction_digest.c dataset.h
	$(CC) $(CFLAGS) -c restriction_digest.c -o restriction_digest_XbaI.o \
 -D_digest_XbaI
restriction_digest_XmnI.o: restriction_digest.c dataset.h
	$(CC) $(CFLAGS) -c restriction_digest.c -o restriction_digest_XmnI.o \
 -D_digest_XmnI
restriction_digest_HindIII.o: restriction_digest.c dataset.h
	$(CC) $(CFLAGS) -c restriction_digest.c \
 -o restriction_digest_HindIII.o -D_digest_HindIII
fasta2kmer: fasta2kmer.c dataset.h kmers.h dataset.o kmers.o binary_array.o
	$(CC) $(CFLAGS) fasta2kmer.c -o ./bin/fasta2kmer dataset.o kmers.o \
 binary_array.o $(PTHREAD) $(MATH)

kmer2pca: kmer2pca.c 
	$(CC) $(CFLAGS) kmer2pca.c -o ./bin/kmer2pca -mavx $(PTHREAD) $(MATH) \
 $(LAPACK)

reverse_with_mask: reverse_with_mask.c dataset.h binary_array.h dataset.o \
                   binary_array.o
	$(CC) $(CFLAGS) reverse_with_mask.c -o ./bin/reverse_with_mask \
 dataset.o binary_array.o $(MATH)

reverse_complement_with_mask: reverse_with_mask.c dataset.h \
                              binary_array.h dataset.o binary_array.o
	$(CC) $(CFLAGS) reverse_complement_with_mask.c -o \
 ./bin/reverse_complement_with_mask dataset.o binary_array.o $(MATH)

digest_XbaI: digester.c restriction_digest.h restriction_digest_XbaI.o \
             dataset.h dataset.o binary_array.h binary_array.o
	$(CC) $(CFLAGS) digester.c -o ./bin/digest_XbaI \
 dataset.o binary_array.o restriction_digest_XbaI.o $(MATH) $(PTHREAD) \
 -D_digest_XbaI

digest_XmnI: digester.c restriction_digest.h restriction_digest_XmnI.o \
             dataset.h dataset.o binary_array.h binary_array.o
	$(CC) $(CFLAGS) digester.c -o ./bin/digest_XmnI \
 dataset.o binary_array.o restriction_digest_XmnI.o $(MATH) $(PTHREAD) \
 -D_digest_XmnI

digest_HindIII: digester.c restriction_digest.h restriction_digest_HindIII.o \
                dataset.h dataset.o binary_array.h binary_array.o
	$(CC) $(CFLAGS) digester.c -o ./bin/digest_HindIII \
 dataset.o binary_array.o restriction_digest_HindIII.o $(MATH) $(PTHREAD) \
 -D_digest_HindIII 

find_satellite: find_satellite.c filter.o filter.h dataset.o binary_array.o \
                dataset.h binary_array.o
	$(CC) $(CFLAGS) find_satellite.c -o ./bin/find_satellite \
 dataset.o binary_array.o filter.o smith_waterman.o $(MATH) $(PTHREAD)

pca2densitymap: pca2densitymap.c dataset.h binary_array.h density.h dataset.o \
                binary_array.o density.o 
	$(CC) $(CFLAGS) pca2densitymap.c -o ./bin/pca2densitymap \
 dataset.o binary_array.o density.o $(PNG) $(MATH)

pca2densityfile: pca2densityfile.c dataset.h binary_array.h density.h dataset.o\
                 binary_array.o density.o cluster.h cluster_io.o volumes.h \
                 volumes.o
	$(CC) $(CFLAGS) pca2densityfile.c -o ./bin/pca2densityfile \
 dataset.o binary_array.o density.o cluster_io.o volumes.o $(PNG) $(MATH)

pca_visual_extract: pca_visual_extract.c dataset.h binary_array.h density.h \
                    cluster.h dataset.o binary_array.o density.o cluster_io.o \
                    volumes.h volumes.o
	$(CC) $(CFLAGS) pca_visual_extract.c -o ./bin/pca_visual_extract \
 dataset.o binary_array.o density.o cluster_io.o volumes.o $(SDL) $(PNG) $(MATH)

virtual_evolution: virtual_evolution.c dataset.h binary_array.h \
                   dataset.o binary_array.o
	$(CC) $(CFLAGS) virtual_evolution.c -o ./bin/virtual_evolution \
 dataset.o binary_array.o $(MATH)

virtual_evolution_controlled: virtual_evolution_controlled.c dataset.h \
                              binary_array.h dataset.o binary_array.o
	$(CC) $(CFLAGS) virtual_evolution_controlled.c -o \
 ./bin/virtual_evolution_controlled dataset.o binary_array.o $(MATH)

simulation_verification: simulation_verification.c dataset.h binary_array.h \
                         cluster.h dataset.o binary_array.o cluster_io.o \
                         volumes.h volumes.o
	$(CC) $(CFLAGS) simulation_verification.c \
 -o ./bin/simulation_verification dataset.o binary_array.o cluster_io.o \
 volumes.o $(MATH)

split_set_from_annotation: split_set_from_annotation.c dataset.h \
                           binary_array.h cluster.h dataset.o binary_array.o \
                           cluster_io.o volumes.o volumes.h
	$(CC) $(CFLAGS) split_set_from_annotation.c \
 -o ./bin/split_set_from_annotation dataset.o binary_array.o cluster_io.o \
 volumes.o $(MATH)

split_set_to_projections: split_set_to_projections.c dataset.h \
                          binary_array.h cluster.h dataset.o binary_array.o \
                          cluster_io.o volumes.h volumes.o 
	$(CC) $(CFLAGS) split_set_to_projections.c \
 -o ./bin/split_set_to_projections dataset.o binary_array.o cluster_io.o \
 volumes.o $(MATH)


find_sequence_in_split_sets: find_sequence_in_split_sets.c dataset.h \
                             binary_array.h cluster.h dataset.o binary_array.o \
                             cluster_io.o volumes.h volumes.o
	$(CC) $(CFLAGS) find_sequence_in_split_sets.c \
 -o ./bin/find_sequence_in_split_sets dataset.o binary_array.o cluster_io.o \
 volumes.o $(MATH)

find_closest_sequence_SW: find_closest_sequence_SW.c dataset.h cluster.h \
                          cluster_io.o dataset.o binary_array.o \
                          smith_waterman.o smith-waterman.h volumes.h volumes.o
	$(CC) $(CFLAGS) find_closest_sequence_SW.c \
 -o ./bin/find_closest_sequence_SW dataset.o binary_array.o cluster_io.o \
 smith_waterman.o volumes.o $(MATH) $(PTHREAD)

tree_map_for_sequence: tree_map_for_sequence.c dataset.h \
                       binary_array.h cluster.h dataset.o binary_array.o \
                       cluster_io.o volumes.o volumes.h 
	$(CC) $(CFLAGS) tree_map_for_sequence.c \
 -o ./bin/tree_map_for_sequence dataset.o binary_array.o cluster_io.o \
 volumes.o $(MATH)

tree_map_for_split_set: tree_map_for_split_set.c dataset.h \
                        binary_array.h cluster.h dataset.o binary_array.o \
                        cluster_io.o volumes.o volumes.h
	$(CC) $(CFLAGS) tree_map_for_split_set.c \
 -o ./bin/tree_map_for_split_set dataset.o binary_array.o cluster_io.o \
 volumes.o $(MATH)


cluster_dbscan_pca: cluster_dbscan_pca.c dbscan.h dataset.h cluster.h \
                    dataset.o cluster_io.o dbscan_L2.o binary_array.o \
                    volumes.h volumes.o
	$(CC) $(CFLAGS) cluster_dbscan_pca.c -o ./bin/cluster_dbscan_pca \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o volumes.o $(MATH) $(PTHREAD)

cluster_dbscan_kmerL1: cluster_dbscan_kmer.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_L1.o binary_array.o \
                       volumes.h volumes.o
	$(CC) $(CFLAGS) cluster_dbscan_kmer.c -o ./bin/cluster_dbscan_kmerL1 \
 dataset.o cluster_io.o dbscan_L1.o binary_array.o volumes.o \
 -D_CLUSTER_KMER_L1 $(MATH) $(PTHREAD)

cluster_dbscan_kmerL2: cluster_dbscan_kmer.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_L2.o binary_array.o \
                       volumes.o volumes.h 
	$(CC) $(CFLAGS) cluster_dbscan_kmer.c -o ./bin/cluster_dbscan_kmerL2 \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o volumes.o \
 -D_CLUSTER_KMER_L2 $(MATH) $(PTHREAD)

cluster_dbscan_SW: cluster_dbscan_SW.c dbscan.h dataset.h cluster.h \
                   dataset.o cluster_io.o dbscan_SW.o binary_array.o \
                   smith_waterman.o volumes.o volumes.h 
	$(CC) $(CFLAGS) cluster_dbscan_SW.c -o ./bin/cluster_dbscan_SW \
 dataset.o cluster_io.o dbscan_SW.o binary_array.o smith_waterman.o volumes.o \
 $(MATH) $(PTHREAD)

cluster_dbscan_SW_GPU: cluster_dbscan_SW.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_SW_GPU.o binary_array.o \
                       volumes.o volumes.h
	$(CC) $(CFLAGS) cluster_dbscan_SW.c -o ./bin/cluster_dbscan_SW_GPU \
 dataset.o cluster_io.o dbscan_SW_GPU.o binary_array.o volumes.o \
 $(MATH) $(OPENCL) -D_SCAN_SMITH_WATERMAN_GPU $(PTHREAD)

adaptive_clustering_SW: adaptive_clustering.c dbscan.h dataset.h cluster.h \
                        dataset.o cluster_io.o dbscan_SW.o binary_array.o \
                        smith_waterman.o volumes.o volumes.h 
	$(CC) $(CFLAGS) adaptive_clustering.c -o ./bin/adaptive_clustering_SW \
 dataset.o cluster_io.o dbscan_SW.o binary_array.o smith_waterman.o volumes.o \
 $(MATH) -D_SCAN_SMITH_WATERMAN $(PTHREAD)

adaptive_clustering_SW_GPU: adaptive_clustering.c dbscan.h dataset.h cluster.h \
	                    dataset.o cluster_io.o dbscan_SW_GPU.o \
                            binary_array.o smith_waterman.o volumes.o volumes.h
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_SW_GPU \
 dataset.o cluster_io.o dbscan_SW_GPU.o binary_array.o smith_waterman.o \
 volumes.o $(MATH) $(OPENCL) -D_SCAN_SMITH_WATERMAN_GPU $(PTHREAD)

adaptive_clustering_SW_MPI_GPU: adaptive_clustering.c dbscan.h dataset.h \
                         cluster.h dataset.o cluster_io.o dbscan_SW_MPI_GPU.o \
                         binary_array.o smith_waterman.o volumes.o volumes.h
	$(MPICC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_SW_MPI_GPU \
 dataset.o cluster_io.o dbscan_SW_MPI_GPU.o binary_array.o smith_waterman.o \
 volumes.o $(MATH) $(OPENCL) $(MPI) -D_SCAN_SMITH_WATERMAN_MPI_GPU $(PTHREAD)

adaptive_clustering_PCA: adaptive_clustering.c dbscan.h dataset.h cluster.h \
	                 dataset.o cluster_io.o dbscan_L2.o binary_array.o \
                         smith_waterman.o volumes.o volumes.h
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_PCA \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o smith_waterman.o \
 volumes.o $(MATH) -D_CLUSTER_PCA $(PTHREAD)

adaptive_clustering_kmer_L2: adaptive_clustering.c dbscan.h dataset.h cluster.h\
	                 dataset.o cluster_io.o dbscan_L2.o binary_array.o \
                         smith_waterman.o volumes.o volumes.h 
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_kmer_L2 \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o smith_waterman.o \
 volumes.o $(MATH) -D_CLUSTER_KMER_L2 $(PTHREAD)

adaptive_clustering_kmer_L1: adaptive_clustering.c dbscan.h dataset.h cluster.h\
	                 dataset.o cluster_io.o dbscan_L1.o binary_array.o \
                         smith_waterman.o volumes.o volumes.h
	$(CC) $(CFLAGS) adaptive_clustering.c -o \
 ./bin/adaptive_clustering_kmer_L1 \
 dataset.o cluster_io.o dbscan_L1.o binary_array.o smith_waterman.o \
 volumes.o $(MATH) -D_CLUSTER_KMER_L1 $(PTHREAD)

split_set_to_fasta: split_set_to_fasta.c dataset.h cluster.h dataset.o \
                    cluster_io.o binary_array.o volumes.o volumes.h
	$(CC) $(CFLAGS) split_set_to_fasta.c -o ./bin/split_set_to_fasta \
 dataset.o cluster_io.o binary_array.o volumes.o $(MATH)

split_sets_to_newick: split_sets_to_newick.c dataset.h cluster.h dataset.o \
                    cluster_io.o binary_array.o volumes.o volumes.h 
	$(CC) $(CFLAGS) split_sets_to_newick.c -o ./bin/split_sets_to_newick \
 dataset.o cluster_io.o binary_array.o volumes.o $(MATH)

filter_split_sets_by_min: filter_split_sets_by_min.c dataset.h cluster.h \
                          dataset.o cluster_io.o binary_array.o volumes.o \
                          volumes.h
	$(CC) $(CFLAGS) filter_split_sets_by_min.c \
 -o ./bin/filter_split_sets_by_min dataset.o cluster_io.o binary_array.o \
 volumes.o $(MATH)

split_set_to_matrix_line: split_set_to_matrix_line.c dataset.h cluster.h \
                          dataset.o cluster_io.o binary_array.o volumes.o \
                          volumes.h
	$(CC) $(CFLAGS) split_set_to_matrix_line.c \
 -o ./bin/split_set_to_matrix_line dataset.o cluster_io.o binary_array.o \
 volumes.o $(MATH)

split_set_to_matrix_annotation: split_set_to_matrix_line.c dataset.h cluster.h \
                                binary_array.o cluster_io.o dataset.o \
                                volumes.o volumes.h
	$(CC) $(CFLAGS) split_set_to_matrix_annotation.c \
 -o ./bin/split_set_to_matrix_annotation dataset.o cluster_io.o \
 binary_array.o volumes.o $(MATH)

print_connections: print_connections.c dataset.h cluster.h dataset.o \
                   cluster_io.o binary_array.o volumes.o volumes.h
	$(CC) $(CFLAGS) print_connections.c -o ./bin/print_connections \
 dataset.o cluster_io.o binary_array.o volumes.o $(MATH)

compareSW: compareSW.c dataset.h comparison.h smith_waterman.o comparison.o \
           dataset.o binary_array.o
	$(CC) $(CFLAGS) compareSW.c -o ./bin/compareSW smith_waterman.o \
 comparison.o dataset.o binary_array.o $(MATH) $(PTHREAD)

silhouette: dataset.h comparison.h dataset.o comparison.o \
            smith_waterman.o binary_array.o silhouette.c
	$(CC) $(CFLAGS) silhouette.c -o ./bin/silhouette smith_waterman.o \
 comparison.o dataset.o binary_array.o $(MATH) $(PTHREAD)


compare_norms: compare_norms.c dataset.h comparison.h dataset.o comparison.o \
               binary_array.o smith_waterman.o
	$(CC) $(CFLAGS) compare_norms.c -o ./bin/compare_norms \
 smith_waterman.o comparison.o dataset.o binary_array.o $(MATH) $(PTHREAD)

consens: consens.c dataset.h dataset.o binary_array.o
	$(CC) $(CFLAGS) consens.c -o ./bin/consens dataset.o binary_array.o \
 $(MATH)

sequence_multiplicity: sequence_multiplicity.c dataset.h binary_array.o \
                       dataset.o
	$(CC) $(CFLAGS) sequence_multiplicity.c -o ./bin/sequence_multiplicity \
 binary_array.o dataset.o $(MATH)

clean:
	rm ./bin/* *.o
