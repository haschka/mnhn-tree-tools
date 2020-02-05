typedef struct {
  void* intensities;
  size_t n_values_per_dimension[3];
  float shift[3];
  float conversion_factor;
  size_t value_size;
  int dimensions;
} density_map;

density_map generate_density_map_from_dataset(dataset ds, int dim,
					      float interval_length);

density_map longmap_to_char_map(density_map map);

int save_2d_density_to_png(char* filename, density_map map);
int save_3d_density_to_pngs(char* filename, density_map map);

void print_2d_density_to_file(FILE* f, density_map map);
