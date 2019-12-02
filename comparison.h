unsigned long* smith_waterman_distances_matrix(dataset ds_one,
					       dataset ds_two);

double mean_from_smith_waterman_distance_matrix(unsigned long* matrix,
						dataset ds_one,
						dataset ds_two);

double sigma_from_smith_waterman_distance_matrix(unsigned long* matix,
						 double mean,
						 dataset ds_one,
						 dataset ds_two);

void shortest_longest_distance_in_matrix(unsigned long* matrix,
					 unsigned long* shortest,
					 unsigned long* largest,
					 dataset ds_one,
					 dataset ds_two);

double silhouette_from_smith_waterman_datasets(dataset ds_one,
					       dataset ds_two);

double sigma_from_smith_waterman_datasets(double mean,
					  dataset ds_one,
					  dataset ds_two);

double mean_from_smith_waterman_datasets(dataset ds_one, dataset ds_two);



