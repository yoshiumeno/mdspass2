#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <fstream>
//#include <iostream>
#include "myheader.h"
#include <cuda.h>
#include <cuda_runtime.h>

__global__ void e_force_morse_cu(int *arrayDint, int *repatom, DOUBLE *rx, DOUBLE *ry, DOUBLE *rz,
				   DOUBLE *fx, DOUBLE *fy, DOUBLE *fz, DOUBLE *epot,
				   DOUBLE *hmat,
				   int *alistnum, int *alist);
void e_force_geam_am();

void potential()
{
  if (strcmp(atom.potential_func, "Morse") == 0) {
    size_t arrayDint_size = sizeof(int) * 2;
    if (istep==0) { cudaMalloc((void **)&arrayDint, arrayDint_size); }
    int arrayHint[2];
    arrayHint[0] = atom.natom; arrayHint[1] = atom.QC;
    cudaMemcpy(arrayDint, arrayHint, arrayDint_size, cudaMemcpyHostToDevice);
    size_t arrayD1_size = sizeof(int) * (atom.natom+1);
    size_t arrayD11_size = sizeof(int) * (NMAX+1);
    size_t arrayD2_size = sizeof(DOUBLE) * (atom.natom+1);
    size_t arrayD3_size = sizeof(DOUBLE) * 9;
    size_t arrayD4_size = sizeof(int) * (NMAX+1)*(NBOOK+1)*4;
    if (istep == 0) {
      cudaMalloc((void **)&arrayDrepatom, arrayD1_size);
      cudaMalloc((void **)&arrayDrx, arrayD2_size);
      cudaMalloc((void **)&arrayDry, arrayD2_size);
      cudaMalloc((void **)&arrayDrz, arrayD2_size);
      cudaMalloc((void **)&arrayDfx, arrayD2_size);
      cudaMalloc((void **)&arrayDfy, arrayD2_size);
      cudaMalloc((void **)&arrayDfz, arrayD2_size);
      cudaMalloc((void **)&arrayDepot, arrayD2_size);
      cudaMalloc((void **)&arrayDhmat, arrayD3_size);
      cudaMalloc((void **)&arrayDalistnum, arrayD11_size);
      cudaMalloc((void **)&arrayDalist, arrayD4_size);
    }
    cudaMemcpy(arrayDrepatom, atom.repatom, arrayD1_size, cudaMemcpyHostToDevice); // arrayD <- arrayH
    cudaMemcpy(arrayDrx, atom.rx, arrayD2_size, cudaMemcpyHostToDevice); // arrayD <- arrayH
    cudaMemcpy(arrayDry, atom.ry, arrayD2_size, cudaMemcpyHostToDevice); // arrayD <- arrayH
    cudaMemcpy(arrayDrz, atom.rz, arrayD2_size, cudaMemcpyHostToDevice); // arrayD <- arrayH
    arrayHhmat[0]=cell.hmat[1][1]; arrayHhmat[1]=cell.hmat[2][1]; arrayHhmat[2]=cell.hmat[3][1];
    arrayHhmat[3]=cell.hmat[1][2]; arrayHhmat[4]=cell.hmat[2][2]; arrayHhmat[5]=cell.hmat[3][2];
    arrayHhmat[6]=cell.hmat[1][3]; arrayHhmat[7]=cell.hmat[2][3]; arrayHhmat[8]=cell.hmat[3][3];
    cudaMemcpy(arrayDhmat, arrayHhmat, arrayD3_size, cudaMemcpyHostToDevice);
    if (istep % book.nbk ==0) {
      cudaMemcpy(arrayDalistnum, book.alistnum, arrayD11_size, cudaMemcpyHostToDevice);
      int icnt=0; //int iarray[(NMAX+1)*(NBOOK+1)*4];
      for (int i=0; i<NMAX+1; i++) {
	for (int j=0; j<NBOOK+1; j++) {
	  for (int k=0; k<4; k++) {
	    iarray[icnt] = book.alist[i][j][k]; icnt++; } } }
      cudaMemcpy(arrayDalist, iarray, arrayD4_size, cudaMemcpyHostToDevice);
    }
    // e_force_morse_cu();
    e_force_morse_cu<<<8, 32>>>(arrayDint, arrayDrepatom, arrayDrx, arrayDry, arrayDrz,
				  arrayDfx, arrayDfy, arrayDfz, arrayDepot,arrayDhmat, arrayDalistnum, arrayDalist);
    cudaMemcpy(atom.fx, arrayDfx, arrayD2_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(atom.fy, arrayDfy, arrayD2_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(atom.fz, arrayDfz, arrayD2_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(atom.epot, arrayDepot, arrayD2_size, cudaMemcpyDeviceToHost);
  } else if (strcmp(atom.potential_func, "GEAM") == 0) {
    e_force_geam_am();
  } else {
    exit(0);
  }
}
