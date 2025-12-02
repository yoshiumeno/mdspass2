#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "myheader.h"
#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif
void bookkeep();
void potential();
void inverse(double a[3][3], double b[3][3]);
void matmul(double a[3][3], double b[3][3], double c[3][3]);
void readconfig(const char* fname);
void readconfig(const char* fname, int mode);
void writeconfig(const char* fname);

extern float neb_tol_fac;
extern int neb_ite, neb_init_read;

void neb_calc(int neb_num, const char* fini, const char* fend)
{
  //int neb_num = 10;
  //int neb_ite = 1000;
  double neb_tol = pow(10,neb_tol_fac);
  double neb_dt = 1e-14;
  double dmax = 0.1*ang;
  double **rxn, **ryn, **rzn, **fxn, **fyn, **fzn, **vxn, **vyn, **vzn;
  double *tx, *ty, *tz, *sx, *sy, *sz, *eneb;
  double qx, qy, qz, xx, yy, zz, dtn;
  double cos0, ff, f0x, f0y, f0z, f1x, f1y, f1z;
  double sp, spr = 5e2, fac = 1e-3;
  double hinmat[3][3];
  dtn = neb_dt;
  writeconfig("CONFIG.TMP");
  readconfig(fini);
  inverse(cell.hmat, hinmat);

  //allocate arrays
  rxn = new double*[atom.natom+1];
  ryn = new double*[atom.natom+1];
  rzn = new double*[atom.natom+1];
  fxn = new double*[atom.natom+1];
  fyn = new double*[atom.natom+1];
  fzn = new double*[atom.natom+1];
  vxn = new double*[atom.natom+1];
  vyn = new double*[atom.natom+1];
  vzn = new double*[atom.natom+1];
  for (int i=0; i<atom.natom+1; i++) {
    rxn[i] = new double[neb_num];
    ryn[i] = new double[neb_num];
    rzn[i] = new double[neb_num];
    fxn[i] = new double[neb_num];
    fyn[i] = new double[neb_num];
    fzn[i] = new double[neb_num];
    vxn[i] = new double[neb_num];
    vyn[i] = new double[neb_num];
    vzn[i] = new double[neb_num];
  }
  tx = new double[atom.natom+1];
  ty = new double[atom.natom+1];
  tz = new double[atom.natom+1];
  sx = new double[atom.natom+1];
  sy = new double[atom.natom+1];
  sz = new double[atom.natom+1];
  eneb = new double[neb_num];

  //reset
  for (int i=0; i<atom.natom+1; i++) {
    for (int j=0; j<neb_num; j++) {
      vxn[i][j] = 0; vyn[i][j] = 0; vzn[i][j] = 0;
      fxn[i][j] = 0; fyn[i][j] = 0; fzn[i][j] = 0;
    }
  }
  double dt2 = dtn/2;
  double dtsq2 = dtn*dt2;

  //set initial nodes
  readconfig(fini,1);
  bookkeep(); potential(); eneb[0] = atom.epotsum;
  for (int i=1; i<=atom.natom; i++) {
    rxn[i][0] = atom.rx[i]; ryn[i][0] = atom.ry[i]; rzn[i][0] = atom.rz[i];
  }
  readconfig(fend,1);
  for (int i=1; i<=atom.natom; i++) {
    rxn[i][neb_num-1] = atom.rx[i]; ryn[i][neb_num-1] = atom.ry[i]; rzn[i][neb_num-1] = atom.rz[i];
  }
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rxn[i][neb_num-1] - rxn[i][0];
    atom.ry[i] = ryn[i][neb_num-1] - ryn[i][0];
    atom.rz[i] = rzn[i][neb_num-1] - rzn[i][0];
    qx = hinmat[0][0]*atom.rx[i] + hinmat[0][1]*atom.ry[i] + hinmat[0][2]*atom.rz[i];
    qy = hinmat[1][0]*atom.rx[i] + hinmat[1][1]*atom.ry[i] + hinmat[1][2]*atom.rz[i];
    qz = hinmat[2][0]*atom.rx[i] + hinmat[2][1]*atom.ry[i] + hinmat[2][2]*atom.rz[i];
    xx = 0; yy = 0; zz = 0;
    if (qx < -0.5) { xx = 1; } if (qx > 0.5) { xx = -1; }
    if (qy < -0.5) { yy = 1; } if (qy > 0.5) { yy = -1; }
    if (qz < -0.5) { zz = 1; } if (qz > 0.5) { zz = -1; }
    rxn[i][neb_num-1] += cell.hmat[0][0]*xx + cell.hmat[0][1]*yy + cell.hmat[0][2]*zz;
    ryn[i][neb_num-1] += cell.hmat[1][0]*xx + cell.hmat[1][1]*yy + cell.hmat[1][2]*zz;
    rzn[i][neb_num-1] += cell.hmat[2][0]*xx + cell.hmat[2][1]*yy + cell.hmat[2][2]*zz;
  }
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rxn[i][neb_num-1]; atom.ry[i] = ryn[i][neb_num-1]; atom.rz[i] = rzn[i][neb_num-1];
  }
  bookkeep(); potential(); eneb[neb_num-1] = atom.epotsum;
  if (neb_init_read == 1) {
    for (int j=1; j<neb_num-1; j++) {
      char filepath[80] = "CONFIG.NEB."; char numc[10];
      snprintf(numc, sizeof(numc), "%02d", j); strcat(filepath,numc);
      readconfig(filepath,1);
      for (int i=1; i<=atom.natom; i++) {
	rxn[i][j] = atom.rx[i]; ryn[i][j] = atom.ry[i]; rzn[i][j] = atom.rz[i];
	atom.rx[i] = rxn[i][j] - rxn[i][0];
	atom.ry[i] = ryn[i][j] - ryn[i][0];
	atom.rz[i] = rzn[i][j] - rzn[i][0];
	qx = hinmat[0][0]*atom.rx[i] + hinmat[0][1]*atom.ry[i] + hinmat[0][2]*atom.rz[i];
	qy = hinmat[1][0]*atom.rx[i] + hinmat[1][1]*atom.ry[i] + hinmat[1][2]*atom.rz[i];
	qz = hinmat[2][0]*atom.rx[i] + hinmat[2][1]*atom.ry[i] + hinmat[2][2]*atom.rz[i];
	xx = 0; yy = 0; zz = 0;
	if (qx < -0.5) { xx = 1; } if (qx > 0.5) { xx = -1; }
	if (qy < -0.5) { yy = 1; } if (qy > 0.5) { yy = -1; }
	if (qz < -0.5) { zz = 1; } if (qz > 0.5) { zz = -1; }
	rxn[i][j] += cell.hmat[0][0]*xx + cell.hmat[0][1]*yy + cell.hmat[0][2]*zz;
	ryn[i][j] += cell.hmat[1][0]*xx + cell.hmat[1][1]*yy + cell.hmat[1][2]*zz;
	rzn[i][j] += cell.hmat[2][0]*xx + cell.hmat[2][1]*yy + cell.hmat[2][2]*zz;
      }
    }
  } else {
    for (int j=1; j<neb_num-1; j++) {
      for (int i=1; i<=atom.natom; i++) {
	xx = (rxn[i][neb_num-1] - rxn[i][0]) / (double)(neb_num - 1);
	yy = (ryn[i][neb_num-1] - ryn[i][0]) / (double)(neb_num - 1);
	zz = (rzn[i][neb_num-1] - rzn[i][0]) / (double)(neb_num - 1);
	rxn[i][j] = rxn[i][0] + xx * (double)j;
	ryn[i][j] = ryn[i][0] + yy * (double)j;
	rzn[i][j] = rzn[i][0] + zz * (double)j;
      }
    }
  }

  std::ofstream feneb("eneb.d");
  std::ofstream feneb_last("eneb_last.d");
  std::ofstream fnebtol("nebtol.d");
  for (int icnt=0; icnt<neb_ite; icnt++) {
    printf("NEB step = %d ", icnt);
    //Check if virtual f and v are not too large
    double scl = 1.0; double scl2 = 1.0;
    xx = 0;
    for (int j=1; j<neb_num-1; j++) {
      for (int i=1; i<=atom.natom; i++) {
	yy = fabs(dtn*vxn[i][j] + dtsq2*fxn[i][j]/atom.wm[i]);
	if (xx < yy) xx = yy;
	yy = fabs(dtn*vyn[i][j] + dtsq2*fyn[i][j]/atom.wm[i]);
	if (xx < yy) xx = yy;
	yy = fabs(dtn*vzn[i][j] + dtsq2*fzn[i][j]/atom.wm[i]);
	if (xx < yy) xx = yy;
      }
    }
    if (xx > dmax) { scl = dmax/xx; scl2 = scl*scl; }
    //vv1
    for (int j=1; j<neb_num-1; j++) {
      for (int i=1; i<=atom.natom; i++) {
	rxn[i][j] += dtn*vxn[i][j]*scl + dtsq2*fxn[i][j]/atom.wm[i]*scl2;
	ryn[i][j] += dtn*vyn[i][j]*scl + dtsq2*fyn[i][j]/atom.wm[i]*scl2;
	rzn[i][j] += dtn*vzn[i][j]*scl + dtsq2*fzn[i][j]/atom.wm[i]*scl2;
	vxn[i][j] += dt2*fxn[i][j]*scl/atom.wm[i];
	vyn[i][j] += dt2*fyn[i][j]*scl/atom.wm[i];
	vzn[i][j] += dt2*fzn[i][j]*scl/atom.wm[i];
      }
    }
    //calculation of force exerted on each node
    for (int j=1; j<neb_num-1; j++) {
      xx = 0;
      for (int i=1; i<=atom.natom; i++) {
	tx[i] = rxn[i][j]-rxn[i][j-1];
	ty[i] = ryn[i][j]-ryn[i][j-1];
	tz[i] = rzn[i][j]-rzn[i][j-1];
	xx += tx[i]*tx[i] + ty[i]*ty[i] + tz[i]*tz[i];
      }
      xx = sqrt(xx);
      yy = 0;
      for (int i=1; i<=atom.natom; i++) {
	tx[i] = rxn[i][j+1]-rxn[i][j];
	ty[i] = ryn[i][j+1]-ryn[i][j];
	tz[i] = rzn[i][j+1]-rzn[i][j];
	yy += tx[i]*tx[i] + ty[i]*ty[i] + tz[i]*tz[i];
      }
      yy = sqrt(yy);
      for (int i=1; i<=atom.natom; i++) {
	tx[i] = (rxn[i][j]-rxn[i][j-1])/xx + (rxn[i][j+1]-rxn[i][j])/yy;
	ty[i] = (ryn[i][j]-ryn[i][j-1])/xx + (ryn[i][j+1]-ryn[i][j])/yy;
	tz[i] = (rzn[i][j]-rzn[i][j-1])/xx + (rzn[i][j+1]-rzn[i][j])/yy;
      }
      zz = 0;
      for (int i=1; i<=atom.natom; i++) {
	zz += (rxn[i][j]-rxn[i][j-1])*(rxn[i][j+1]-rxn[i][j])
	  +   (ryn[i][j]-ryn[i][j-1])*(ryn[i][j+1]-ryn[i][j])
	  +   (rzn[i][j]-rzn[i][j-1])*(rzn[i][j+1]-rzn[i][j]);
      }
      cos0 = zz/xx/yy;
      //normalize tangent vector
      //xx = 0; // <-----?
      for (int i=1; i<=atom.natom; i++) {
	xx += tx[i]*tx[i] + ty[i]*ty[i] + tz[i]*tz[i];
      }
      xx = sqrt(xx);
      for (int i=1; i<=atom.natom; i++) {
	tx[i] /= xx; ty[i] /= xx; tz[i] /= xx;
      }
      //calc dV (gradient of potential)
      for (int i=1; i<=atom.natom; i++) {
	atom.rx[i] = rxn[i][j]; atom.ry[i] = ryn[i][j]; atom.rz[i] = rzn[i][j];
      }
      potential(); eneb[j] = atom.epotsum;
      //calc Fs (spring force)
      for (int i=1; i<=atom.natom; i++) {
	sx[i] = spr * (rxn[i][j+1] - rxn[i][j] * 2 + rxn[i][j-1]);
	sy[i] = spr * (ryn[i][j+1] - ryn[i][j] * 2 + ryn[i][j-1]);
	sz[i] = spr * (rzn[i][j+1] - rzn[i][j] * 2 + rzn[i][j-1]);
      }
      //calc force
      xx = 0; yy = 0;
      for (int i=1; i<=atom.natom; i++) {
	xx += atom.fx[i]*tx[i] + atom.fy[i]*ty[i] + atom.fz[i]*tz[i];
	yy +=      sx[i]*tx[i] +      sy[i]*ty[i] +      sz[i]*tz[i];
      }
      ff = 0.5 * (1.0 + cos(M_PI * cos0));
      for (int i=1; i<=atom.natom; i++) {
	f0x = atom.fx[i] - xx*tx[i] + yy*tx[i];
	f0y = atom.fy[i] - xx*ty[i] + yy*ty[i];
	f0z = atom.fz[i] - xx*tz[i] + yy*tz[i];
	f1x = f0x + ff * (sx[i] - yy*tx[i]);
	f1y = f0y + ff * (sy[i] - yy*ty[i]);
	f1z = f0z + ff * (sz[i] - yy*tz[i]);
	fxn[i][j] = 0; fyn[i][j] = 0; fzn[i][j] = 0;
	if (atom.mfx[i] == false) { fxn[i][j] = f1x; }
	if (atom.mfy[i] == false) { fyn[i][j] = f1y; }
	if (atom.mfz[i] == false) { fzn[i][j] = f1z; }
      }
    } // j = 1 ~ j<neb_num-1
    for (int j=0; j<neb_num; j++) {
      feneb << j << " " << eneb[j]/ev << std::endl;
    }

    //optimization
    //vv2
    for (int j=1; j<neb_num-1; j++) {
      for (int i=1; i<=atom.natom; i++) {
	vxn[i][j] += dt2*fxn[i][j]/atom.wm[i];
	vyn[i][j] += dt2*fyn[i][j]/atom.wm[i];
	vzn[i][j] += dt2*fzn[i][j]/atom.wm[i];
      }
    }
    //damper
    for (int j=1; j<neb_num-1; j++) {
      sp = 0;
      for (int i=1; i<=atom.natom; i++) {
	sp += fxn[i][j]*vxn[i][j]+fyn[i][j]*vyn[i][j]+fzn[i][j]*vzn[i][j];
      }
      if (sp < 0 ) {
	for (int i=1; i<=atom.natom; i++) {
	  vxn[i][j] = 0; vyn[i][j] = 0; vzn[i][j] = 0;
	}
      }
    }
    //convergence check
    xx = 0;
    for (int j=1; j<neb_num-1; j++) {
      for (int i=1; i<=atom.natom; i++) {
	yy = fxn[i][j]*fxn[i][j]; if (yy > xx) { xx = yy; }
	yy = fyn[i][j]*fyn[i][j]; if (yy > xx) { xx = yy; }
	yy = fzn[i][j]*fzn[i][j]; if (yy > xx) { xx = yy; }
      }
    }
    fnebtol << icnt << " " << sqrt(xx) << std::endl;
    printf("Residual error = %e  Tolerance = %e\n",sqrt(xx),neb_tol);
    if (sqrt(xx) < neb_tol) { printf("NEB done\n"); break; }
    //printf("%e is still larger than %e\n",sqrt(xx),neb_tol);
    feneb << std::endl;
  } // loop icnt end

  { using namespace std;
    for (int j=0; j<neb_num; j++) {
      feneb_last <<setw(4)<< j << " " << setiosflags(ios::fixed) << setw(10) 
		 <<setprecision(6)<< eneb[j]/ev << endl;
    }
  }
  feneb.close(); fnebtol.close(); feneb_last.close();

  for (int j=1; j<neb_num-1; j++) {
    for (int i=1; i<=atom.natom; i++) {
      atom.rx[i] = rxn[i][j]; atom.ry[i] = ryn[i][j]; atom.rz[i] = rzn[i][j];
    }
    char filepath[80] = "CONFIG.NEB."; char numc[10];
    snprintf(numc, sizeof(numc), "%02d", j); strcat(filepath,numc);
    writeconfig(filepath);
  }

  //deallocate arrays
  for (int i=0; i<atom.natom+1; i++) {
    delete[] rxn[i], ryn[i], rzn[i];
    delete[] fxn[i], fyn[i], fzn[i];
    delete[] vxn[i], vyn[i], vzn[i];
  }
  delete[] rxn, ryn, rzn, fxn, fyn, fzn, vxn, vyn, vzn;
  delete[] tx, ty, tz, sx, sy, sz, eneb;

  readconfig("CONFIG.TMP");
}
