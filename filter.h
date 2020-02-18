dataset filter_ds_by_size(dataset in_set, size_t size,
			  size_t plus, size_t minus);


dataset filter_by_SW_distance(dataset in_set,
			      char* target,
			      size_t target_length,
			      unsigned int max_distance,
			      int n_threads);
