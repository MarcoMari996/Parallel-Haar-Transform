#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define RED "\e[0;31m"

/***
 * definition of stuctures useful for the computations
 */
typedef unsigned char uint8_t;

struct haar_coefficients {
    uint8_t *approximation;
    uint8_t *horizontal_details;
    uint8_t *vertical_details;
    uint8_t *diagonal_details;
    int height;
    int width;
};

struct decomposition {
    int levels;
    struct haar_coefficients *coefficients;
};

unsigned char * load_image_from_text(char* filename, int *height, int *width) {
    FILE *image_file = fopen(filename, "r");
    if(image_file == NULL)
    {
        fprintf(stderr, RED "[ERROR] there's a problem with the image file!\n\n");
        exit(8);
    }
    int h, w;
    fscanf(image_file, "%d\t%d", &h, &w);
    *height = h;
    *width = w;
    unsigned char *image = (unsigned char*)malloc(h*w*sizeof(unsigned char));
    int i, j;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            fscanf(image_file, "%u\t", &image[i * w + j]);
        }
    }

    fclose(image_file);
    return image;
}

uint8_t * inverse_haar_transform(struct haar_coefficients coef, int height, int width) {
    uint8_t *result = (uint8_t*) malloc(height*width*sizeof(uint8_t));
    for (int i=0; i < height/2; i++) {
        for (int j=0; j < width/2; j++) {
            result[2*i*width + 2*j] = coef.approximation[i*width/2 + j] + coef.horizontal_details[i*width/2 + j] +
                                      coef.vertical_details[i*width/2 + j] + coef.diagonal_details[i*width/2 + j];
            result[2*i*width + 2*j+1] = coef.approximation[i*width/2 + j] - coef.horizontal_details[i*width/2 + j] +
                                        coef.vertical_details[i*width/2 + j] - coef.diagonal_details[i*width/2 + j];
            result[(2*i+1)*width + 2*j] = coef.approximation[i*width/2 + j] + coef.horizontal_details[i*width/2 + j] -
                                          coef.vertical_details[i*width/2 + j] - coef.diagonal_details[i*width/2 + j];
            result[(2*i+1)*width + 2*j+1] = coef.approximation[i*width/2 + j] - coef.horizontal_details[i*width/2 + j] -
                                            coef.vertical_details[i*width/2 + j] + coef.diagonal_details[i*width/2 + j];
        }
    }
    return result;
}

struct haar_coefficients haar_transform(uint8_t *img, int height, int width) {
    struct haar_coefficients coefficients;
    coefficients.approximation = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.horizontal_details = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.vertical_details = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.diagonal_details = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.height = height/2;
    coefficients.width = width/2;

    int i, j, tid;
    int chunk = coefficients.height / omp_get_num_threads(); // set loop iteration chunk size //
    // printf("\nChunk size: %d\n", chunk);

    #pragma omp parallel private(tid, i,j) shared(height, width, img, chunk)
    {
        #pragma omp for schedule(auto)
           for (i = 0; i < height; i += 2) {
               for (j = 0; j < width; j += 2) {
                   // tid = omp_get_thread_num();
                   // printf("Thread %d is doing:\ti = %d\tj = %d\n", tid, i, j);
                   uint8_t a = img[i * width + j];
                   uint8_t b = img[i * width + (j + 1)];
                   uint8_t c = img[(i + 1) * width + j];
                   uint8_t d = img[(i + 1) * width + (j + 1)];
                   // low pass on rows - low pass on cols --> approximation
                   coefficients.approximation[i / 2 * width / 2 + j / 2] = (uint8_t) (((a + b) + (c + d)) / 4);
                   // low pass on rows - high pass on cols --> horizontal details
                   coefficients.horizontal_details[i / 2 * width / 2 + j / 2] = (uint8_t) (((a + b) - (c + d)) / 4);
                   // high pass on rows - low pass on cols --> vertical details
                   coefficients.vertical_details[i / 2 * width / 2 + j / 2] = (uint8_t) (((a - b) + (c - d)) / 4);
                   // high pass on rows - high pass on cols --> diagonal details
                   coefficients.diagonal_details[i / 2 * width / 2 + j / 2] = (uint8_t) (((a - b) - (c - d)) / 4);
               }
           }
    } // end of the paralell region

    return coefficients;
}

struct decomposition compute_decomposition(uint8_t *img, int height, int width, int levels) {
    struct decomposition decomp;
    decomp.levels = levels;
    decomp.coefficients = (struct haar_coefficients*) malloc(decomp.levels * sizeof(struct haar_coefficients) );

    decomp.coefficients[0] = haar_transform(img, height, width);
    for (int l=1; l < levels; l++) {
        decomp.coefficients[l] = haar_transform(decomp.coefficients[l - 1].approximation, decomp.coefficients[l - 1].height, decomp.coefficients[l - 1].width);
    }
    return decomp;
}

void check_applicability(unsigned char * img, int height, int width, int levels) {
    if (img == NULL) {
        printf(RED "Error in loading the image\n");
        exit(1);
    }
    if ((height % 2 != 0) | (width % 2 != 0)) {
        printf(RED "Load an image whose sizes are a multiple of 2 !");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    double begin = omp_get_wtime();
    if (argc != 3) {
        fprintf(stderr, RED "[ERROR]: usage of the command \n\t haar [levels] [image file path]\n\n");
        exit(9);
    }
    int levels = atoi(argv[1]);
    char *filename = argv[2];
    int height, width;
    unsigned char *img = load_image_from_text(filename, &height,  &width);
    check_applicability(img, height, width, levels);

    // compute decomposition levels
    omp_set_num_threads(24);
    struct decomposition decomposition;
    decomposition = compute_decomposition(img, height, width, levels);

    // image reconstruction
    // uint8_t *reconstructed_img = inverse_haar_transform(haar_coefficients, height, width);

    // print wall-clock time
    double end = omp_get_wtime();
    printf("Exec time: \t%.3f ms\n", (end - begin)*1000);
    return 0;
}

