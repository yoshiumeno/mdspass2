#include <GL/glui.h>
#include <png.h>
#include <string>
#include <stdlib.h>
//#include <iostream>
//#include <fstream>
//#include <math.h>

#if !defined __linux__ && !defined __APPLE__
#define snprintf    _snprintf
#pragma comment(lib, "libpng15.lib")
#endif

extern int capture_count;
void capture(const char* fname);

void capture()
{
  capture("___");
}

//void capture()
void capture(const char* fname)
{
  bool fname_given = true;
  if (strcmp(fname, "___")==0 ) { fname_given = false; }
#ifndef NOPNG
  char filepath[80] = "SNAPSHOT";
  if (!fname_given) {
    if (capture_count<0) { capture_count = 0; }
    int num = capture_count;
    char numc[10];
    if (num<10000) { snprintf(numc, sizeof(numc), "%04d", num);
    } else { snprintf(numc, sizeof(numc), "%d", num); }
    strcat(filepath,numc); strcat(filepath,".png");
  } else {
    strcpy(filepath,fname);
  }
    png_bytep raw1D;
    png_bytepp raw2D;
    int i;
    int width = glutGet(GLUT_WINDOW_WIDTH);
    int height = glutGet(GLUT_WINDOW_HEIGHT);

    // 構造体確保
    FILE *fp = fopen(filepath, "wb");
    png_structp pp = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop ip = png_create_info_struct(pp);
    // 書き込み準備
    png_init_io(pp, fp);
    png_set_IHDR(pp, ip, width, height,
        8, // 8bit以外にするなら変える
        PNG_COLOR_TYPE_RGBA, // RGBA以外にするなら変える
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    // ピクセル領域確保
    raw1D = (png_bytep)malloc(height * png_get_rowbytes(pp, ip));
    raw2D = (png_bytepp)malloc(height * sizeof(png_bytep));
    for (i = 0; i < height; i++)
        raw2D[i] = &raw1D[i*png_get_rowbytes(pp, ip)];
    // 画像のキャプチャ
#if defined __linux__ || defined __APPLE__
    glReadBuffer(GL_FRONT);
#else
	glReadBuffer(GL_BACK);
#endif
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 初期値は4
    glReadPixels(0, 0, width, height,
            GL_RGBA, // RGBA以外にするなら変える
            GL_UNSIGNED_BYTE, // 8bit以外にするなら変える
            (void*)raw1D);
    // 上下反転
    for (i = 0; i < height/ 2; i++){
        png_bytep swp = raw2D[i];
        raw2D[i] = raw2D[height - i - 1];
        raw2D[height - i - 1] = swp;
    }
    // 書き込み
    png_write_info(pp, ip);
    png_write_image(pp, raw2D);
    png_write_end(pp, ip);
    // 開放
    png_destroy_write_struct(&pp, &ip);
    fclose(fp);
    free(raw1D);
    free(raw2D);

    printf("Screen is captured to '%s'\n", filepath);
    if (!fname_given) { capture_count++; }
    //glui->sync_live();
#endif
}
