//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "myheader.h"
#include <cuda.h>
#include <cuda_runtime.h>
//DOUBLE v(DOUBLE rmeter);
//DOUBLE vp(DOUBLE rmeter);

//__global__ void e_force_morse_cuda(DOUBLE *array, int len)
__global__ void e_force_morse_cu(int *arrayDint, int *repatom, DOUBLE *rx, DOUBLE *ry, DOUBLE *rz,
				   DOUBLE *fx, DOUBLE *fy, DOUBLE *fz, DOUBLE *epot,
				   DOUBLE *hmat,
				   int *alistnum, int *ialist)
{
  DOUBLE rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  int natom = arrayDint[0]; int QC = arrayDint[1];
  DOUBLE hmat1[4], hmat2[4], hmat3[4];
  hmat1[1]=hmat[0]; hmat2[1]=hmat[1]; hmat3[1]=hmat[2];
  hmat1[2]=hmat[3]; hmat2[2]=hmat[4]; hmat3[2]=hmat[5];
  hmat1[3]=hmat[6]; hmat2[3]=hmat[7]; hmat3[3]=hmat[8];

  for (int i=1; i<=natom; i++) {
    fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;  epot[i] = 0.0;
  }

  DOUBLE ep = 0.3429; DOUBLE al = 1.3588; DOUBLE ro = 2.866; //DOUBLE ev = 1.6021892e-19;
  DOUBLE rcut = 13.0e-10; DOUBLE rcut2 = rcut * rcut;
  int tid = threadIdx.x + blockDim.x * blockIdx.x;
  int nthreads = blockDim.x * gridDim.x;
  int part = natom / nthreads; if (natom%nthreads!=0) { part++; }
  int imin = part*tid+1; int imax = part*(tid+1); if (imax>natom) {imax=natom;}
  //for (int i=1; i<=natom; i++) {
  for (int i=imin; i<=imax; i++) {
    if (alistnum[i]>0) {
      for (int k=1; k<=alistnum[i]; k++) {
	//j  = alist[i][k][0];
	j = ialist[(NBOOK+1)*4*i+4*k+0];
	//if (j>=i) {
	if ((j>=i)||(j<i)) { // <== Double counting algo
	  if ((QC==0)||(repatom[i]==1)||(repatom[j]==1)) {
	    //ix = alist[i][k][1]; iy = alist[i][k][2]; iz = alist[i][k][3];
	    ix = ialist[(NBOOK+1)*4*i+4*k+1];
	    iy = ialist[(NBOOK+1)*4*i+4*k+2];
	    iz = ialist[(NBOOK+1)*4*i+4*k+3];
	    // rr2 = Dist2(i,j,ix,iy,iz);
	    DOUBLE xx = rx[j] + hmat1[1]*ix+hmat1[2]*iy+hmat1[3]*iz - rx[i];
	    DOUBLE yy = ry[j] + hmat2[1]*ix+hmat2[2]*iy+hmat2[3]*iz - ry[i];
	    DOUBLE zz = rz[j] + hmat3[1]*ix+hmat3[2]*iy+hmat3[3]*iz - rz[i];
	    rr2 = ( xx*xx + yy*yy + zz*zz );
	    if (rr2 < rcut2) {
	      rr = sqrt(rr2);
	      // drx = Dx(i,j,ix,iy,iz); dry = Dy(i,j,ix,iy,iz); drz = Dz(i,j,ix,iy,iz);
	      drx = rx[j] + hmat1[1]*ix+hmat1[2]*iy+hmat1[3]*iz - rx[i];
	      dry = ry[j] + hmat2[1]*ix+hmat2[2]*iy+hmat2[3]*iz - ry[i];
	      drz = rz[j] + hmat3[1]*ix+hmat3[2]*iy+hmat3[3]*iz - rz[i];
	      // vp0 = vp(rr); v0 = v(rr);
	      DOUBLE rang = rr * 1.0e10;
	      vp0 = -2.0*al*ep*(exp(-2.0*al*(rang-ro))-exp(-al*(rang-ro))) *ev*1.0e10;
	      v0 = ep*(exp(-2.0*al*(rang-ro))-2.0*exp(-al*(rang-ro))) *ev;

	      if ((QC==0)||(repatom[i]==1)) {
		fx[i] = fx[i]+vp0/rr*drx;
		fy[i] = fy[i]+vp0/rr*dry;
		fz[i] = fz[i]+vp0/rr*drz;
		epot[i]=epot[i]+v0/2.0;
		//if (i==j) { epot[i]=epot[i]+v0/2.0; } // <== Double counting algo
	      }
	      if ((QC==0)||(repatom[j]==1)) {
		//fx[j] = fx[j]-vp0/rr*drx; // <== Double counting algo
		//fy[j] = fy[j]-vp0/rr*dry;
		//fz[j] = fz[j]-vp0/rr*drz;
		//epot[j]=epot[j]+v0/2.0;
	      }

	    }
	  }
	} //j>=i
      }
    }
    //    printf("FAST %d %20.15e %20.15e\n",i,fx[i],rr2);
  }
  //  if (istep==0) {exit(0);}

}

