#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include <SDL2/SDL.h>                                                           

#include"dataset.h"
#include"density.h"
#include"cluster.h"

#define SDL_MAIN_HANDLED

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}

union pixel {
  unsigned int i;
  unsigned char c[4];
};

int pos_check_start_x_lt_endx_start_y_lt_endy(int pos_x, int pos_y,
					      int start_x, int start_y,
					      int end_x, int end_y) {
  return(start_x <= pos_x && pos_x <= end_x &&
	 start_y <= pos_y && pos_y <= end_y);
}

int pos_check_start_x_gt_endx_start_y_lt_endy(int pos_x, int pos_y,
					      int start_x, int start_y,
					      int end_x, int end_y) {
  return(start_x >= pos_x && pos_x >= end_x &&
	 start_y <= pos_y && pos_y <= end_y);
}

int pos_check_start_x_lt_endx_start_y_gt_endy(int pos_x, int pos_y,
					      int start_x, int start_y,
					      int end_x, int end_y) {
  return(start_x <= pos_x && pos_x <= end_x &&
	 start_y >= pos_y && pos_y >= end_y);
}

int pos_check_start_x_gt_endx_start_y_gt_endy(int pos_x, int pos_y,
					      int start_x, int start_y,
					      int end_x, int end_y) {
  return(start_x >= pos_x && pos_x >= end_x &&
	 start_y >= pos_y && pos_y >= end_y);
}


cluster get_cluster_from_rectangle(dataset ds, int dim_a,int dim_b,
				   float min_a,float  min_b,
				   float  max_a, float max_b,
				   int width, int height,
				   int rect_start_x,int rect_start_y,
				   int x, int y, float delta) {

  int i;
  int pos_x,pos_y;

  int (*checker) (int, int, int, int, int, int);

  cluster cl;
    
  cl.members = (int*)malloc(1000*sizeof(int));
  cl.n_members = 0;
  
  if(rect_start_x <= x && rect_start_y <= y) {
    checker = &pos_check_start_x_lt_endx_start_y_lt_endy;
  } else if (rect_start_x >= x && rect_start_y <= y) {
    checker = &pos_check_start_x_gt_endx_start_y_lt_endy;
  } else if (rect_start_x <= x && rect_start_y >= y) {
    checker = &pos_check_start_x_lt_endx_start_y_gt_endy;
  } else if (rect_start_x >= x && rect_start_y >= y) {
    checker = &pos_check_start_x_gt_endx_start_y_gt_endy;
  }
  
  for(i = 0; i<ds.n_values; i++) {
    pos_x = (int)fabsf((ds.values[dim_a][i]-min_a)/delta);
    pos_y = (int)fabsf((ds.values[dim_b][i]-min_b)/delta);
    if(checker(pos_x,pos_y,rect_start_x,rect_start_y,x,y)) {
      cl.n_members++;
      if(cl.n_members % 1000 == 0) {
	cl.members = (int*)realloc(cl.members,sizeof(int)*(cl.n_members+1000));
      }
      cl.members[cl.n_members-1] = i;
    }
  }
  if(cl.n_members != 0) {
    cl.members = (int*)realloc(cl.members,sizeof(int)*cl.n_members);
  } else {
    free(cl.members);
  }
  return(cl);
}

static inline int position_in_image(int width,
				    int height,
				    int x, int y) {
  return((x > 0 && x < width && y > 0 && y < height));
}

static inline void draw_rect(unsigned int* frame, int start_x, int start_y,
			     int end_x, int end_y, int width, int height) {


  int i;
  int s_x, e_x;
  int s_y, e_y;

  union pixel red;
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  red.c[0] = 0;
  red.c[1] = 0;
  red.c[2] = 255;
  red.c[3] = 0;
#else
  red.c[3] = 0;
  red.c[2] = 0;
  red.c[1] = 255;
  red.c[0] = 0;
#endif
  
  if(start_x > end_x) {
    s_x = end_x; e_x = start_x;
  } else {
    s_x = start_x; e_x = end_x;
  }

  if(start_y > end_y) {
    s_y = end_y; e_y = start_y;
  } else {
    s_y = start_y; e_y = end_y;
  }

  for(i = s_x; i<e_x;i++) {
    frame[width*s_y+i] = red.i;
    frame[width*e_y+i] = red.i;
  }
  for(i = s_y;i<e_y;i++) {
    frame[width*i+s_x] = red.i;
    frame[width*i+e_x] = red.i;
  }
} 
    

int main(int argc, char** argv) {

  dataset ds;

  int i;

  size_t input_dimensions;

  FILE* fasta_f;
  FILE* projection_f;

  float* min_max;

  float a_min, a_max, b_min, b_max, len_a, len_b, ratio, delta;

  int width, height;
  int pos_x,pos_y;

  unsigned int* image_frame;
  unsigned int* static_pca_image;

  union pixel lit;

  SDL_Window *image_window = NULL;
  SDL_Renderer *image_renderer = NULL;
  SDL_Texture *image_texture;

  SDL_Event event;
  
  int pitch;

  int mouse_down = 0;
  int save_avail = 0;
  
  int rect_start_x;
  int rect_start_y;

  cluster cl;

  char* save_file_name = argv[6];

  int cycles = 0;
  int dim_a, dim_b;

  int x,y;

  if(argc < 6) {
    printf("Arguments are \n");
    printf(" [file-in] fasta file to select sequences from\n");
    printf(" [file-in] file containing pca projections corresponding to the \n"
	   "           fasta file \n");
    printf(" (int)     dimensions in pca projections file \n");
    printf(" (int)     first dimension for 2D representation \n");
    printf(" (int)     second dimension for 2D representation \n");
    printf(" [file-out] fasta file that will contain the selected sequences\n");
    return(1);
  }

  sscanf(argv[3],"%lu",&input_dimensions);
  sscanf(argv[4],"%i",&dim_a);
  sscanf(argv[5],"%i",&dim_b);

  dim_a--;
  dim_b--;
  
  cl.members = NULL;
  cl.n_members = 0;

  
  lit.c[0] = 255;
  lit.c[1] = 255;
  lit.c[2] = 255;
  lit.c[3] = 255;
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (projection_f = fopen(argv[2], "r"))) file_error(argv[2]);
  load_projections_from_file_into_dataset(projection_f,input_dimensions,&ds);
  fclose(projection_f);

  min_max = get_min_max_in_dimension_from_dataset(ds, dim_a);
  a_min = min_max[0];
  a_max = min_max[1];
  free(min_max);

  min_max = get_min_max_in_dimension_from_dataset(ds, dim_b);
  b_min = min_max[0];
  b_max = min_max[1];
  free(min_max);
  
  len_a = fabsf(a_max-a_min);
  len_b = fabsf(b_max-b_min);

  delta = len_a/399.f;

  width = 400;
  height = (int)(len_b/delta+1.f);

  pitch = width*4;
  
  static_pca_image = (unsigned int*)malloc(sizeof(unsigned int)*width*height);

  memset(static_pca_image,0,sizeof(unsigned int)*width*height);

  for(i = 0 ; i < ds.n_values; i++ ) {
    pos_x = (int)(fabsf(ds.values[dim_a][i]-a_min)/delta);
    pos_y = (int)(fabsf(ds.values[dim_b][i]-b_min)/delta);

    static_pca_image[pos_y*width+pos_x] = lit.i;

  }

  SDL_SetMainReady();
  SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS);

  image_window = SDL_CreateWindow("Fancy PCA Sequence Selector!",
				  100, 20, width, height, 0);
  
  
#ifdef SOFTWARE_RENDERING
  image_renderer = SDL_CreateRenderer(image_window,-1,SDL_RENDERER_SOFTWARE);
#else
  image_renderer = SDL_CreateRenderer(image_window,
				      -1,SDL_RENDERER_ACCELERATED);
#endif

  image_texture = SDL_CreateTexture(image_renderer, SDL_PIXELFORMAT_ARGB8888,
				    SDL_TEXTUREACCESS_STREAMING,
				    width,height);

  SDL_LockTexture(image_texture,NULL,(void**)&image_frame, &pitch);

  memcpy(image_frame,static_pca_image,
	 sizeof(unsigned int)*width*height);

  SDL_UnlockTexture(image_texture);
  
  SDL_RenderCopy(image_renderer,image_texture,NULL,NULL);
  SDL_RenderPresent(image_renderer);


  while(1) {

    while(SDL_PollEvent(&event)) {
      switch(event.type) {
	
      case SDL_QUIT:
	goto finish;
	break;

      case SDL_MOUSEBUTTONDOWN:
	rect_start_x = event.button.x;
	rect_start_y = event.button.y;
	mouse_down = 1;

	break;

      case SDL_MOUSEBUTTONUP:
	printf("Calculating... \n");
	if(cl.n_members != 0) {
	  free(cl.members);
	  cl.members = 0;
	}
	cl = get_cluster_from_rectangle(ds, dim_a, dim_b,
					a_min, b_min, a_max, b_max,
					width, height,
					rect_start_x, rect_start_y,
					event.button.x, event.button.y, delta);
	if(cl.n_members > 0) {
	  printf("Selected %i Sequences ! \n",cl.n_members);
	  printf("Hit enter to save to fasta file: %s \n",save_file_name);
	  save_avail = 1;
	} else {
	  printf("No Sequences in selection - Try again! \n",cl.n_members);
	}
	mouse_down = 0;
	
	break;

      case SDL_KEYDOWN:
	if(event.key.keysym.sym == SDLK_RETURN ||
	   event.key.keysym.sym == SDLK_RETURN2) {
	  if(save_avail) {
	    printf("Saving...\n");
	    create_single_cluster_file(save_file_name, cl,ds);
	    printf("Saved !\n");
	    save_avail = 0;
	  }
	}
	break;
	
      case SDL_WINDOWEVENT:
	if(event.window.event == SDL_WINDOWEVENT_CLOSE) {
	  goto finish;
	}
	break;
      }
    }

    if (mouse_down) {
      /*
      x = event.motion.x;
      y = event.motion.y;
      */
      SDL_GetMouseState(&x,&y);
      if(position_in_image(width,height,x,y)) {

	SDL_LockTexture(image_texture,NULL,(void**)&image_frame, &pitch);

	memcpy(image_frame,static_pca_image,
	       sizeof(unsigned int)*width*height);

	draw_rect(image_frame,rect_start_x,rect_start_y, x,y,width, height);

	SDL_UnlockTexture(image_texture);
	SDL_RenderCopy(image_renderer,image_texture,NULL,NULL);
	
	SDL_RenderPresent(image_renderer);

      }
    }

    if ( cycles == 128 ) {
      SDL_RenderPresent(image_renderer);
      cycles = 0;
    }

    SDL_Delay(10);
    cycles++;
  }
 finish:
  
  free(static_pca_image);
  if(cl.n_members != 0) free(cl.members);
  free_values_from_dataset(ds);
  free_sequences_from_dataset(ds);
  
  SDL_DestroyTexture(image_texture);
  SDL_DestroyRenderer(image_renderer);
  SDL_DestroyWindow(image_window);

  SDL_Quit();
  return(0);
}
