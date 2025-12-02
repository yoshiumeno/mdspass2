#include <iostream>
#include <stdio.h>
#include <fstream>
#include "myheader.h"
void potential();
void inverse(double a[3][3], double b[3][3]);
void matmul(double a[3][3], double b[3][3], double c[3][3]);
void multiply_cell(int ix, int iy, int iz, int mode);
void divide_cell(int ix, int iy, int iz, int mode);
void forceconst(double **fcmat);
void kgen_fcc(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_bcc(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_hcp(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_sc(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_chain(double *kx, double *ky, double *kz, int phonon_knum, int &size);

struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;
extern "C"
void zgeev_( const char* jobvl, const char* jobvr, int* n, dcomplex* a,
	     int* lda, dcomplex* w, dcomplex* vl, int* ldvl, dcomplex* vr, int* ldvr,
	     dcomplex* work, int* lwork, double* rwork, int* info);

extern int phonon_rep, phonon_knum, phonon_kp;

void phonon_calc()
{
  //printf("Phonon calculation\n");
  //multiply cell size
  int n0 = atom.natom*3;
  int ix = phonon_rep, iy = phonon_rep, iz = phonon_rep;
  multiply_cell(ix, iy, iz, 1);
  //allocate array
  double **fcmat;
  dcomplex *dm;
  int nn = atom.natom*3;
  fcmat = new double*[nn];
  dm = new dcomplex[n0*n0];
  for (int i=0; i<nn; i++) {
    fcmat[i] = new double[nn];
  }
  int nn2 = n0*2, info = 1, ldvl = 1, ldvr = 1;
  dcomplex *work; work = new dcomplex[nn2];
  double *rwork; rwork = new double[nn2];
  dcomplex *w;
  dcomplex *vl, *vr;
  w = new dcomplex[n0];
  vl = new dcomplex[n0*ldvl]; vr = new dcomplex[nn*ldvr];

  forceconst(fcmat);
  std::ofstream ofs("phonon.d");
  printf("Phonon dispersion data is written in phonon.d\n");

  double *kx, *ky, *kz;

  int size = 0;
  if (phonon_kp == 0) {
    kgen_fcc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 1) {
    kgen_bcc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 2) {
    kgen_sc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 3) {
    kgen_hcp(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 4) {
    kgen_chain(kx, ky, kz, phonon_knum, size);
  }
  printf("# of k-points = %d\n",size);
  kx = new double[size];
  ky = new double[size];
  kz = new double[size];
  if (phonon_kp == 0) {
    kgen_fcc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 1) {
    kgen_bcc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 2) {
    kgen_sc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 3) {
    kgen_hcp(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 4) {
    kgen_chain(kx, ky, kz, phonon_knum, size);
  }

  for (int ik=0; ik<size; ik++) {

    for (int i=0; i<n0*n0; i++) {
      dm[i].re = 0; dm[i].im = 0; }
    
    for (int i=0; i<n0; i++) {
      int ia = i/3 + 1;
      for (int j=0; j<nn; j++) {
	int ja = j/3 + 1;
	double dx = atom.rx[ia] - atom.rx[ja];
	double dy = atom.ry[ia] - atom.ry[ja];
	double dz = atom.rz[ia] - atom.rz[ja];
	double xx = - (kx[ik]*dx + ky[ik]*dy + kz[ik]*dz);
	double ms = sqrt(atom.wm[ia]*atom.wm[ja]);
	double xxr = fcmat[i][j] * cos(xx) / ms;
	double xxi = fcmat[i][j] * sin(xx) / ms;
	int ii = i+(j%n0)*n0;
	dm[ii].re += xxr; dm[ii].im += xxi;
      }
    }
    zgeev_("N", "N", &n0, dm, &n0, w, vl, &ldvl, vr, &ldvr, work, &nn2, rwork, &info);
    
    double scale = 1e-27;
    for (int i=0; i<n0; i++) {
      ofs << ik << " " << w[i].re * scale << std::endl;
    }
    ofs << std::endl;
  }

  //deallocate array
  for (int i=0; i<nn; i++) {
    delete[] fcmat[i];
  }
  delete[] fcmat, dm, work, rwork, w, vl, vr;
  delete[] kx, ky, kz;
  ofs.close();

  //Divide cell
  divide_cell(ix, iy, iz, 1);

}

void kgen_fcc(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[5] = {1,1,1,0.5,1}; // ratio of k-point density
  double spx[10],spy[10],spz[10];
  double a;
  if (cell.hmat[1][0]/cell.hmat[0][0] > 0.5) {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  } else {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  }
  const char *sp[] = {"W", "L", "L", "Gamma", "Gamma", "X", "X", "U", "K", "Gamma"};
  printf("lattice constant a = %f (A)\n",a/ang);
  spx[0] = 1.00; spy[0] = 0.00; spz[0] = 0.50; // W
  spx[1] = 0.50; spy[1] = 0.50; spz[1] = 0.50; // L
  spx[2] = 0.50; spy[2] = 0.50; spz[2] = 0.50; // L
  spx[3] = 0.00; spy[3] = 0.00; spz[3] = 0.00; // Gamma
  spx[4] = 0.00; spy[4] = 0.00; spz[4] = 0.00; // Gamma
  spx[5] = 1.00; spy[5] = 0.00; spz[5] = 0.00; // X
  spx[6] = 1.00; spy[6] = 0.00; spz[6] = 0.00; // X
  spx[7] = 1.00; spy[7] = 0.25; spz[7] = 0.25; // U
  spx[8] = 0.75; spy[8] = 0.00; spz[8] = 0.75; // K
  spx[9] = 0.00; spy[9] = 0.00; spz[9] = 0.00; // Gamma
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<10; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[9]; ky[size]=spy[9]; kz[size]=spz[9];
  }
  size++;
}

void kgen_bcc(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[5] = {1,1,0.5,1,1}; // ratio of k-point density
  double spx[10],spy[10],spz[10];
  double a;
  if (fabs(cell.hmat[0][1]/cell.hmat[0][0]) > 0.5) {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  } else {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  }
  const char *sp[] = {"Gamma", "H", "H", "N", "N", "P", "P", "Gamma", "Gamma", "N"};
  printf("lattice constant a = %e\n",a);
  spx[0] = 0.00; spy[0] = 0.00; spz[0] = 0.00; // Gamma
  spx[1] = 1.00; spy[1] = 0.00; spz[1] = 0.00; // H
  spx[2] = 1.00; spy[2] = 0.00; spz[2] = 0.00; // H
  spx[3] = 0.50; spy[3] = 0.00; spz[3] = 0.50; // N
  spx[4] = 0.50; spy[4] = 0.00; spz[4] = 0.50; // N
  spx[5] = 0.50; spy[5] = 0.50; spz[5] = 0.50; // P
  spx[6] = 0.50; spy[6] = 0.50; spz[6] = 0.50; // P
  spx[7] = 0.00; spy[7] = 0.00; spz[7] = 0.00; // Gamma
  spx[8] = 0.00; spy[8] = 0.00; spz[8] = 0.00; // Gamma
  spx[9] = 0.50; spy[9] = 0.00; spz[9] = 0.50; // N
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<10; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[9]; ky[size]=spy[9]; kz[size]=spz[9];
  }
  size++;
}

void kgen_sc(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[4] = {1,1,1,1}; // ratio of k-point density
  double spx[8],spy[8],spz[8];
  double a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  const char *sp[] = {"Gamma", "X", "X", "M", "M", "R", "R", "Gamma"};
  printf("lattice constant a = %e\n",a);
  spx[0] = 0.00; spy[0] = 0.00; spz[0] = 0.00; // Gamma
  spx[1] = 1.00; spy[1] = 0.00; spz[1] = 0.00; // X
  spx[2] = 1.00; spy[2] = 0.00; spz[2] = 0.00; // X
  spx[3] = 1.00; spy[3] = 1.00; spz[3] = 0.00; // M
  spx[4] = 1.00; spy[4] = 1.00; spz[4] = 0.00; // M
  spx[5] = 1.00; spy[5] = 1.00; spz[5] = 1.00; // R
  spx[6] = 1.00; spy[6] = 1.00; spz[6] = 1.00; // R
  spx[7] = 0.00; spy[7] = 0.00; spz[7] = 0.00; // Gamma
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<10; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[7]; ky[size]=spy[7]; kz[size]=spz[7];
  }
  size++;
}

void kgen_hcp(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[7] = {1,0.5,1,0.5,1,0.5,1}; // ratio of k-point density
  double spx[14],spy[14],spz[14];
  double a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  double c = cell.hmat[2][2]/(double)(phonon_rep*2-1);
  const char *sp[] = {"Gamma", "M", "M", "K", "K", "Gamma", "Gamma", "A", "A", "L", "L", "H", "H", "A"};
  printf("lattice constant a = %e\n",a);
  spx[ 0] = 0.00; spy[ 0] = 0.00; spz[ 0] = 0.00; // Gamma
  spx[ 1] = 1.00; spy[ 1] = 0.00; spz[ 1] = 0.00; // M
  spx[ 2] = 1.00; spy[ 2] = 0.00; spz[ 2] = 0.00; // M
  spx[ 3] = 1.00; spy[ 3] = 0.50; spz[ 3] = 0.00; // K
  spx[ 4] = 1.00; spy[ 4] = 0.50; spz[ 4] = 0.00; // K
  spx[ 5] = 0.00; spy[ 5] = 0.00; spz[ 5] = 0.00; // Gamma
  spx[ 6] = 0.00; spy[ 6] = 0.00; spz[ 6] = 0.00; // Gamma
  spx[ 7] = 0.00; spy[ 7] = 0.00; spz[ 7] = 1.00; // A
  spx[ 8] = 0.00; spy[ 8] = 0.00; spz[ 8] = 1.00; // A
  spx[ 9] = 1.00; spy[ 9] = 0.00; spz[ 9] = 1.00; // L
  spx[10] = 1.00; spy[10] = 0.00; spz[10] = 1.00; // L
  spx[11] = 1.00; spy[11] = 0.50; spz[11] = 1.00; // H
  spx[12] = 1.00; spy[12] = 0.50; spz[12] = 1.00; // H
  spx[13] = 0.00; spy[13] = 0.00; spz[13] = 1.00; // A
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<14; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/c* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[13]; ky[size]=spy[13]; kz[size]=spz[13];
  }
  size++;
}

void kgen_chain(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[1] = {1}; // ratio of k-point density
  double spx[2],spy[2],spz[2];
  double a;
  a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  printf("lattice constant a = %e\n",a);
  spx[0] = 0.00; spy[0] = 0.00; spz[0] = 0.00; // Gamma
  spx[1] = 1.00; spy[1] = 0.00; spz[1] = 0.00; // X
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<2; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[1]; ky[size]=spy[1]; kz[size]=spz[1];
  }
  size++;
}
