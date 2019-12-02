CC=gcc
CFLAGS=-O2 -march=native
LAPACK=-llapack
MATH=-lm
PTHREAD=-pthread

all: fasta2kmer kmer2pca cluster_dbscan_pca cluster_dbscan_kmerL1 \
     cluster_dbscan_kmerL2 cluster_dbscan_SW compareSW silhouette

compare.o: compare.c compare.h dataset.h smith-waterman.h
	$(CC) $(CFLAGS) -c comparison.c -o comparison
smith_waterman.o: smith-waterman.c smith-waterman.h
	$(CC) $(CFLAGS) -c smith-waterman.c -o smith_waterman.o
binary_array.o: binary_array.c binary_array.h
	$(CC) $(CFLAGS) -c binary_array.c -o binary_array.o
cluster_io.o: cluster_io.c cluster.h dataset.h
	$(CC) $(CFLAGS) -c cluster_io.c -o cluster_io.o
dbscan_SW.o: dbscan.c dbscan.h cluster.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_SMITH_WATERMAN -o dbscan_SW.o
dbscan_L1.o: dbscan.c dbscan.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_L1 -o dbscan_L1.o
dbscan_L2.o: dbscan.c dbscan.h cluster.h
	$(CC) $(CFLAGS) -c dbscan.c -D_SCAN_L2 -o dbscan_L2.o 
dataset.o: dataset.c dataset.h
	$(CC) $(CFLAGS) -c dataset.c -o dataset.o
kmers.o: kmers.c dataset.h kmers.h
	$(CC) $(CFLAGS) -c kmers.c -o kmers.o

fasta2kmer: fasta2kmer.c dataset.h kmers.h dataset.o kmers.o
	$(CC) $(CFLAGS) fasta2kmer.c -o ./bin/fasta2kmer dataset.o kmers.o \
 $(PTHREAD)

kmer2pca: kmer2pca.c 
	$(CC) $(CFLAGS) kmer2pca.c -o ./bin/kmer2pca -mavx $(PTHREAD) $(MATH) \
 $(LAPACK)

cluster_dbscan_pca: cluster_dbscan_pca.c dbscan.h dataset.h cluster.h \
                    dataset.o cluster_io.o dbscan_L2.o binary_array.o
	$(CC) $(CFLAGS) cluster_dbscan_pca.c -o ./bin/cluster_dbscan_pca \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o

cluster_dbscan_kmerL1: cluster_dbscan_kmer.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_L1.o binary_array.o
	$(CC) $(CFLAGS) cluster_dbscan_kmer.c -o ./bin/cluster_dbscan_kmerL1 \
 dataset.o cluster_io.o dbscan_L1.o binary_array.o -D_CLUSTER_KMER_L1

cluster_dbscan_kmerL2: cluster_dbscan_kmer.c dbscan.h dataset.h cluster.h \
                       dataset.o cluster_io.o dbscan_L2.o binary_array.o
	$(CC) $(CFLAGS) cluster_dbscan_kmer.c -o ./bin/cluster_dbscan_kmerL2 \
 dataset.o cluster_io.o dbscan_L2.o binary_array.o -D_CLUSTER_KMER_L2

cluster_dbscan_SW: cluster_dbscan_SW.c dbscan.h dataset.h cluster.h \
                   dataset.o cluster_io.o dbscan_SW.o binary_array.o \
                   smith_waterman.o
	$(CC) $(CFLAGS) cluster_dbscan_SW.c -o ./bin/cluster_dbscan_SW \
 dataset.o cluster_io.o dbscan_SW.o binary_array.o smith_waterman.o

compareSW: compareSW.c dataset.h comparison.h smith_waterman.o comparison.o \
           dataset.o
	$(CC) $(CFLAGS) compareSW.c -o ./bin/compareSW smith_waterman.o \
 comparison.o dataset.o $(MATH)

silhouette: dataset.h comparison.h dataset.h dataset.o comparison.o \
            smith_waterman.o 
	$(CC) $(CFLAGS) silhouette.c -o ./bin/silhouette smith_waterman.o \
 comparison.o dataset.o $(MATH)

clean:
	rm ./bin/* *.o
