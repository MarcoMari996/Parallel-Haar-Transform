#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

/***
 * definition of stuctures useful for the computations
 */
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
        printf("Error with the image file!");
        exit(1);
    }
    int h, w;
    fscanf(image_file, "%d\t%d\n", &h, &w);
    unsigned char *image = (unsigned char*)malloc(h*w*sizeof(unsigned char));
    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++) {
            fscanf(image_file, "%u\t", &image[i * w + j]);
        }
    }

    *height = h;
    *width = w;

    fclose(image_file);
    return image;
}

/***
 * Converting the original image into gray levels
 * @param img the original image
 * @param height
 * @param width
 * @param channels
 * @return the image in gray levels
 */
uint8_t * to_gray_levels(unsigned char *img, int height, int width, int channels) {
    size_t img_size = width * height * channels;
    int gray_channels = channels == 4 ? 2 : 1;
    size_t gray_img_size = width * height * gray_channels;
    unsigned char *gray_img = malloc(gray_img_size);
    if (gray_img == NULL) {
        printf("unable to allocate memory for the image. \n");
        exit(1);
    }
    for(unsigned char *p = img, *pg = gray_img; p != img + img_size; p += channels, pg += gray_channels) {
        *pg = (uint8_t) ((*p + *(p + 1) + *(p + 2)) / 3.0);
        if (channels == 4) {
            *(pg + 1) = *(p + 3);
        }
    }
    return gray_img;
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

uint8_t * merge_coefficients(struct haar_coefficients coef, int height, int width) {
    uint8_t *res = (uint8_t*) malloc(width*height*sizeof(uint8_t));
    int i,j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            if (i < height/2) {
                if (j < width / 2) {
                    res[i * width + j] = coef.approximation[i * width/2 + j];
                } else {
                    res[i * width + j] = coef.vertical_details[i * width/2 + j-height/2];
                }
            } else {
                if (j < width / 2) {
                    res[i * width + j] = coef.horizontal_details[(i-height/2) * width/2 + j];
                } else {
                    res[i * width + j] = coef.diagonal_details[(i - height/2) * width/2 + j-height/2];
                }
            }
        }
    }
    return res;
}

uint8_t * merge_decomposition(struct decomposition d) {
    uint8_t *res;
    for (int l=d.levels-1; l>0; l--) {
        d.coefficients[l-1].approximation = merge_coefficients(d.coefficients[l], 2*d.coefficients[l].height, 2*d.coefficients[l].width);
    }
    res = merge_coefficients(d.coefficients[0], 2*d.coefficients[0].height, 2*d.coefficients[0].width);
    return res;
}

struct haar_coefficients haar_transform(uint8_t *img, int height, int width) {
    struct haar_coefficients coefficients;
    coefficients.approximation = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.horizontal_details = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.vertical_details = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.diagonal_details = (uint8_t *) malloc(height / 2 * width / 2 * sizeof(uint8_t));
    coefficients.height = height/2;
    coefficients.width = width/2;
    for (int i = 0; i < height; i += 2) {
        for (int j = 0; j < width; j += 2) { // applying the filter to all the columns
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

void check_applicability(unsigned char * img, int height, int width) {
    if (img == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    if ((height % 2 != 0) | (width % 2 != 0)) {
        printf("Load an image whose sizes are a multiple of 2 !");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    clock_t begin = clock();
    // // should be provided by terminal
    char *filename = "/images/horse2048.txt";
    int height, width;
    unsigned char *img = load_image_from_text(filename, &height,  &width);
    // // should be provided by terminal
    int levels = 3;
    // int width, height, channels;
    // // loading the image
    // unsigned char *img = stbi_load(filename,&width, &height, &channels, 0);
    // // limiting applicability on sizes multiple of 2
    check_applicability(img, height, width);
    // // convert image into gray levels
    // uint8_t *gray_img = to_gray_levels(img, height, width, channels);

    // compute decomposition levels
    struct decomposition decomposition;

    decomposition = compute_decomposition(img, height, width, levels);

    // // produce final output
    // uint8_t *merged_details = merge_decomposition(decomposition);
    // stbi_write_jpg("/out.jpeg", width, height, 1, merged_details, 100);
    // // trhesholding ??
    // // uint8_t *reconstructed_img = inverse_haar_transform(haar_coefficients, height, width);

    // // freeing memory allocations
    // stbi_image_free(img);
    // stbi_image_free(gray_img);

    // print wall-clock time
    clock_t end = clock();
    printf("Exec time: \t%.3f ms\n", (double)(end - begin) * 1000 / CLOCKS_PER_SEC);

    return 0;
}

