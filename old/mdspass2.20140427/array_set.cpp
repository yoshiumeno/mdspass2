#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"
#include <GL/glui.h>

extern GLfloat **color, **color0;
extern int *iatom, *repidx;

void deallocate_arrays()
{
 // Deallocate arrays
  if (color) { for (int i=1; i<=atom.natom*3; i++) { if (color[i]) free(color[i]); } free(color); }
  if (color0) { for (int i=1; i<=atom.natom*3; i++) { if (color0[i]) free(color0[i]); } free(color0); }
  if (iatom) { free(iatom); }
  if (repidx) { free(repidx); }
  if (atom.asp) { for (int i=1; i<=atom.natom; i++) { if (atom.asp[i]) free(atom.asp[i]); } free(atom.asp);  }
  if (atom.anum) {delete[] atom.anum; atom.anum = NULL; }
  if (atom.wm) { delete[] atom.wm; atom.wm = NULL;}
  if (atom.rx) {free(atom.rx);}; if (atom.ry) {free(atom.ry);}; if (atom.rz) {free(atom.rz);}
  if (atom.rx_float) {free(atom.rx_float);}; if (atom.ry_float) {free(atom.ry_float);}; if (atom.rz_float) {free(atom.rz_float);}
  if (atom.rx_org) {free(atom.rx_org);}; if (atom.ry_org) {free(atom.ry_org);}; if (atom.rz_org) {free(atom.rz_org);}
  if (atom.rx_p) {free(atom.rx_p);}; if (atom.ry_p) {free(atom.ry_p);}; if (atom.rz_p) {free(atom.rz_p);}
  if (atom.vx) {free(atom.vx);}; if (atom.vy) {free(atom.vy);}; if (atom.vz) {free(atom.vz);}
  if (atom.fx) {free(atom.fx);}; if (atom.fy) {free(atom.fy);}; if (atom.fz) {free(atom.fz);}
  if (atom.ax) { delete[] atom.ax; atom.ax = NULL; }
  if (atom.ay) { delete[] atom.ay; atom.ay = NULL; }
  if (atom.az) { delete[] atom.az; atom.az = NULL; }
  if (atom.bx) { delete[] atom.bx; atom.bx = NULL; }
  if (atom.by) { delete[] atom.by; atom.by = NULL; }
  if (atom.bz) { delete[] atom.bz; atom.bz = NULL; }
  if (atom.cx) { delete[] atom.cx; atom.cx = NULL; }
  if (atom.cy) { delete[] atom.cy; atom.cy = NULL; }
  if (atom.cz) { delete[] atom.cz; atom.cz = NULL; }
  if (atom.fx_l) { delete[] atom.fx_l; atom.fx_l = NULL; }
  if (atom.fy_l) { delete[] atom.fy_l; atom.fy_l = NULL; }
  if (atom.fz_l) { delete[] atom.fz_l; atom.fz_l = NULL; }
  if (atom.fx_float) {free(atom.fx_float);}; if (atom.fy_float) {free(atom.fy_float);}; if (atom.fz_float) {free(atom.fz_float);}
  if (atom.mfx) {free(atom.mfx);}; if (atom.mfy) {free(atom.mfy);}; if (atom.mfz) {free(atom.mfz);}
  if (atom.lock) {free(atom.lock);};
  if (atom.epot) {free(atom.epot);}
  if (atom.epot_p) {free(atom.epot_p);}
  if (atom.epot_float) {free(atom.epot_float);}
  if (atom.repatom) {free(atom.repatom);}
  if (atom.elem_id) {free(atom.elem_id);}
  if (atom.evecx) { for (int i=1; i<=atom.natom; i++) { free(atom.evecx[i]); } free(atom.evecx);}
  if (atom.evecy) { for (int i=1; i<=atom.natom; i++) { free(atom.evecy[i]); } free(atom.evecy);}
  if (atom.evecz) { for (int i=1; i<=atom.natom; i++) { free(atom.evecz[i]); } free(atom.evecz);}
  if (atom.satom) {
    for (int i=1; i<=atom.natom; i++) {
      for (int j=0; j<3; j++) {	free(atom.satom[i][j]); }  free(atom.satom[i]); } 
    free(atom.satom); }
  if (atom.nneighbor) { delete[] atom.nneighbor; atom.nneighbor = NULL;
    for (int i=0; i<=atom.natom; i++) { delete[] atom.neighbor[i]; }
    delete[] atom.neighbor; atom.neighbor = NULL; }
  delete[] book.alistnum;
  if (cg.sgx) { for (int i=0; i<=atom.natom; i++) { delete[] cg.sgx[i]; }
    delete[] cg.sgx; cg.sgx = NULL; }
  if (cg.sgy) { for (int i=0; i<=atom.natom; i++) { delete[] cg.sgy[i]; }
    delete[] cg.sgy; cg.sgy = NULL; }
  if (cg.sgz) { for (int i=0; i<=atom.natom; i++) { delete[] cg.sgz[i]; }
    delete[] cg.sgz; cg.sgz = NULL; }
  if (cg.shx) { for (int i=0; i<=atom.natom; i++) { delete[] cg.shx[i]; }
    delete[] cg.shx; cg.shx = NULL; }
  if (cg.shy) { for (int i=0; i<=atom.natom; i++) { delete[] cg.shy[i]; }
    delete[] cg.shy; cg.shy = NULL; }
  if (cg.shz) { for (int i=0; i<=atom.natom; i++) { delete[] cg.shz[i]; }
    delete[] cg.shz; cg.shz = NULL; }
  if (cg.rstx) { delete[] cg.rstx; cg.rstx = NULL; }
  if (cg.rsty) { delete[] cg.rsty; cg.rsty = NULL; }
  if (cg.rstz) { delete[] cg.rstz; cg.rstz = NULL; }
 // Deallocate arrays end
}

void allocate_arrays()
{
  // Allocate arrays
  //  color = (GLfloat **)malloc(sizeof(GLfloat*)*(atom.natom*3+1));
  //  for (int i=1; i<=atom.natom*3; i++) { color[i] = (GLfloat *)malloc(sizeof(GLfloat)*4); }
  color = (GLfloat **)malloc(sizeof(GLfloat*)*(atom.natom*3+1));
  color0 = (GLfloat **)malloc(sizeof(GLfloat*)*(atom.natom*3+1));
  for (int i=1; i<=atom.natom*3; i++) { color[i] = (GLfloat *)malloc(sizeof(GLfloat)*4); }
  for (int i=1; i<=atom.natom*3; i++) { color0[i] = (GLfloat *)malloc(sizeof(GLfloat)*4); }
  iatom  = (int *)calloc(atom.natom*2+1, sizeof(int));
  repidx = (int *)calloc(atom.natom*2+1, sizeof(int));
  atom.asp = (char **)malloc(sizeof(char*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.asp[i] = (char *)malloc(sizeof(char)*3); }
  atom.anum = new int[atom.natom+1];
  atom.wm = new double[atom.natom+1];
  atom.rx = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ry = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rz = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ax = new double[atom.natom+1];
  atom.ay = new double[atom.natom+1];
  atom.az = new double[atom.natom+1];
  atom.bx = new double[atom.natom+1];
  atom.by = new double[atom.natom+1];
  atom.bz = new double[atom.natom+1];
  atom.cx = new double[atom.natom+1];
  atom.cy = new double[atom.natom+1];
  atom.cz = new double[atom.natom+1];
  atom.rx_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.ry_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.rz_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.rx_org = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ry_org = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rz_org = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rx_p = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ry_p = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rz_p = (double *)calloc(atom.natom+1, sizeof(double));
  atom.vx = (double *)calloc(atom.natom+1, sizeof(double));
  atom.vy = (double *)calloc(atom.natom+1, sizeof(double));
  atom.vz = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fx = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fy = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fz = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fx_l = new double[atom.natom+1];
  atom.fy_l = new double[atom.natom+1];
  atom.fz_l = new double[atom.natom+1];
  atom.fx_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.fy_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.fz_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.mfx = (bool *)calloc(atom.natom+1, sizeof(bool));
  atom.mfy = (bool *)calloc(atom.natom+1, sizeof(bool));
  atom.mfz = (bool *)calloc(atom.natom+1, sizeof(bool));
  atom.lock= (bool *)calloc(atom.natom+1, sizeof(bool));
  atom.epot = (double *)calloc(atom.natom+1, sizeof(double));
  atom.epot_p = (double *)calloc(atom.natom+1, sizeof(double));
  atom.epot_float = (FLOAT2 *)calloc(atom.natom+1, sizeof(FLOAT2));
  atom.repatom = (int *)calloc(atom.natom+1, sizeof(int));
  atom.elem_id = (int *)calloc(atom.natom+1, sizeof(int));
  atom.evecx = (double **)malloc(sizeof(double*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.evecx[i] = (double *)malloc(sizeof(double)*MAXMODE); }
  atom.evecy = (double **)malloc(sizeof(double*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.evecy[i] = (double *)malloc(sizeof(double)*MAXMODE); }
  atom.evecz = (double **)malloc(sizeof(double*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.evecz[i] = (double *)malloc(sizeof(double)*MAXMODE); }
  atom.satom = (double ***)malloc(sizeof(double**)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) {
    atom.satom[i] = (double **)malloc(sizeof(double*)*3);
    for (int j=0; j<3; j++) { atom.satom[i][j] = (double *)malloc(sizeof(double)*3); } }
  atom.nneighbor = new int[atom.natom+1]; atom.neighbor = new int*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) { atom.neighbor[i] = new int[MAXNEIGHBOR]; }
  book.alistnum = new int[atom.natom+1]; book.alloc = true;
  
  cg.sgx = new double*[atom.natom+1];
  cg.sgy = new double*[atom.natom+1];
  cg.sgz = new double*[atom.natom+1];
  cg.shx = new double*[atom.natom+1];
  cg.shy = new double*[atom.natom+1];
  cg.shz = new double*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) {
    cg.sgx[i] = new double[2];
    cg.sgy[i] = new double[2];
    cg.sgz[i] = new double[2];
    cg.shx[i] = new double[2];
    cg.shy[i] = new double[2];
    cg.shz[i] = new double[2];
  }
  cg.rstx = new double[atom.natom+1];
  cg.rsty = new double[atom.natom+1];
  cg.rstz = new double[atom.natom+1];
  // Allocate arrays end
}
