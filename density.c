#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<unistd.h>

#include<float.h>
#include<png.h>

#include"dataset.h"
#include"density.h"

static inline float* get_min_max_in_dimension_from_dataset(dataset ds,
							   int dim) {
  int i;

  float min = FLT_MAX;
  float max = FLT_MIN;

  float * v = ds.values[dim];

  float * min_max = (float*)malloc(sizeof(float)*2);
  
  for( i = 0 ; i < ds.n_values ; i++ ) {
    if( v[i] < min ) { min = v[i]; }
    if( v[i] > max ) { max = v[i]; }
  }

  min_max[0] = min;
  min_max[1] = max;
  return(min_max);
}

density_map generate_density_map_from_dataset(dataset ds, int dim,
					      float interval_length) {

  int i,j;

  float* min_max[3];

  density_map map;
  
  size_t n_intensities;
  unsigned long* intensities;

  size_t offset, dimension_multiplier;
  
  size_t w_index, h_index, z_index;
  
  map.dimensions = dim;

  n_intensities = 1;
  for (i=0; i<dim; i++) {
    min_max[i] = get_min_max_in_dimension_from_dataset(ds,i);
    map.n_values_per_dimension[i] =
      (size_t)((float)fabsf(min_max[i][1]-min_max[i][0])
	       /(float)interval_length)+1;
    n_intensities *= map.n_values_per_dimension[i];
    map.shift[i] = 0 - min_max[i][0];
  }
  
  intensities = (unsigned long*)malloc(sizeof(unsigned long)*n_intensities);
  memset(intensities,0,sizeof(unsigned long)*n_intensities);

  
  
  for (i=0; i<ds.n_values; i++) {
    offset = 0;
    dimension_multiplier = 1;
    switch (dim) {
    case 2:
      w_index = (size_t)((fabsf((float)ds.values[0][i]
				+(float)map.shift[0])/(float)interval_length));
      h_index = (size_t)((fabsf((float)ds.values[1][i]
				+(float)map.shift[1])/(float)interval_length));
      intensities[h_index*map.n_values_per_dimension[0]
		  +w_index]++;
      break;
    case 3:
      w_index = (size_t)((fabsf((float)ds.values[0][i]
				+(float)map.shift[0])/(float)interval_length));
      h_index = (size_t)((fabsf((float)ds.values[1][i]
				+(float)map.shift[1])/(float)interval_length));
      z_index = (size_t)((fabsf((float)ds.values[2][i]
				+(float)map.shift[2])/(float)interval_length));
      intensities[z_index
		  *map.n_values_per_dimension[1]*map.n_values_per_dimension[0]
		  +h_index
		  *map.n_values_per_dimension[0]
		  +w_index]++;
      break;
    default:
      printf("Dimension error!\n");
      _exit(1);
    }
  }
  
  map.conversion_factor = 1;
  map.value_size = sizeof(unsigned long);
  map.intensities = (void*) intensities;
  return(map);
}

density_map longmap_to_char_map(density_map map) {

  int i;
  size_t n_intensities;
  
  unsigned long max;
  
  density_map cmap;

  unsigned char* c_intensities;
  unsigned long* l_intensities = (unsigned long*)map.intensities;

  float factor;
  
  n_intensities = 1;
  for (i=0; i<map.dimensions; i++) {
    n_intensities*= map.n_values_per_dimension[i];
  }

  c_intensities = (unsigned char*)malloc(sizeof(unsigned char)*n_intensities);

  max = 0;

  for (i=0; i<n_intensities;i++) {
    if ( max < l_intensities[i]) {
      max = l_intensities[i];
    }
  }

  factor = (float)255.f/(float)max;
  
  for(i=0; i<n_intensities;i++) {
    /*
    if (l_intensities[i] > 0) {
      c_intensities[i]=255;
    } else {
      c_intensities[i]=0;
    };
    */
    if ( factor*l_intensities[i] > 255 ) {
      c_intensities[i] = 255 ;
    } else {
      c_intensities[i] = factor*l_intensities[i];
    }
  }

  cmap.intensities = (void*)c_intensities;

  for(i=0; i<map.dimensions; i++) {
    cmap.n_values_per_dimension[i] = map.n_values_per_dimension[i];
    cmap.shift[i] = map.shift[i];
  }
  cmap.conversion_factor = factor;
  cmap.value_size = sizeof(unsigned char);
  cmap.dimensions = map.dimensions;

  return(cmap);
}

void print_2d_density_to_file(FILE* f, density_map map) {
  size_t i,j;

  unsigned long* intensities = (unsigned long*)map.intensities;

  if (map.value_size != sizeof(unsigned long)) {
    printf("2d density print only works for \"unsigned long\" density maps");
    _exit(1);
  }

  for(j = 0; j < map.n_values_per_dimension[1]; j++) {
    for(i = 0; i < map.n_values_per_dimension[0]; i++) {

      fprintf(f,"%lu %lu %lu\n",
	      i,j,intensities[j*map.n_values_per_dimension[0]+i]);
    }
  }
}


int save_2d_density_to_png(char* filename, density_map map) {

  int i,j;
  
  png_bytep *row = NULL;
  png_bytep png_pixel;

  FILE* f = fopen(filename,"wb");

  png_structp structure = png_create_write_struct(PNG_LIBPNG_VER_STRING,
						  NULL, NULL, NULL);
 
  unsigned char* pixel_c = (unsigned char*)map.intensities;

  png_infop info;

  unsigned int width = map.n_values_per_dimension[0];
  unsigned int height = map.n_values_per_dimension[1];

  if (f == NULL ) return 1;
  if (structure == NULL) return 2;

  info = png_create_info_struct(structure);

  png_init_io(structure,f);
  
  png_set_IHDR(structure, info, width, height, 8,
	       PNG_COLOR_TYPE_GRAY,
	       PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);
  
  png_write_info(structure, info);
  
  row = (png_bytep*)malloc(sizeof(png_bytep)*height);
  
  if(row == NULL) return 3;

  for(i = 0; i < height; i++) {
    row[i] = (png_byte*)malloc(sizeof(unsigned char)*width);
    if(row[i] == 0) return 3; 
  }

  for(j = 0; j < height; j++) {
    for(i = 0; i < width; i++) {
      png_pixel = (row[j]+i);
      png_pixel[0] = pixel_c[j*width+i];
    }
  }

  png_write_image(structure, row);
  png_write_end(structure, NULL);

  for(i = 0; i < height; i++) {
    free(row[i]);
  }
  free(row);

  png_destroy_write_struct(&structure, &info);
  
  fclose(f);
}

int save_3d_density_to_pngs(char* filename, density_map map) {

  int i,j,k;
  
  png_bytep *row = NULL;
  png_bytep png_pixel;

  FILE* f;

  png_structp structure;

  unsigned char* pixel_c = (unsigned char*)map.intensities;

  png_infop info;

  unsigned int width = map.n_values_per_dimension[0];
  unsigned int height = map.n_values_per_dimension[1];

  unsigned int stack = map.n_values_per_dimension[3];

  size_t offset;
  
  char final_file_name[256];
  char char_k[5];

  if (k > 9999) {
    return(4);
  }
  
  for (k = 0 ; k < stack ; k++) {

    for(i=0;i<5;i++) {
      char_k[i] = 0;
    }

    sprintf(char_k,"%04d",k);
    final_file_name[0] = 0;
    strcat(final_file_name,filename);
    strcat(final_file_name,char_k);
    strcat(final_file_name,".png");
    
    f = fopen(final_file_name,"wb");

    structure = png_create_write_struct(PNG_LIBPNG_VER_STRING,
					NULL, NULL, NULL);

    offset = k*width*height;

    if (f == NULL ) return 1;
    if (structure == NULL) return 2;
    
    png_init_io(structure,f);
  
    png_set_IHDR(structure, info, width, height, 8,
		 PNG_COLOR_TYPE_GRAY,
		 PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT,
		 PNG_FILTER_TYPE_DEFAULT);

    png_write_info(structure, info);

    row = (png_bytep*)malloc(sizeof(png_bytep)*height);

    if(row == NULL) return 3;

    for(i = 0; i < height; i++) {
      row[i] = (png_byte*)malloc(sizeof(unsigned char)*width);
      if(row[i] == 0) return 3; 
    }

    for(j = 0; j < height; j++) {
      for(i = 0; i < width; i++) {
	png_pixel = (row[j]+i);
	png_pixel[0] = pixel_c[offset+j*width+i];
      }
    }

    png_write_image(structure, row);
    png_write_end(structure, NULL);

    for(i = 0; i < height; i++) {
      free(row[i]);
    }
    free(row);

    png_destroy_write_struct(&structure, &info);
     
    fclose(f);
  }
}
