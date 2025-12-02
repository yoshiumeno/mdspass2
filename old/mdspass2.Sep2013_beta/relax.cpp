#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

void vscale();
void vscale(double target);
void vscale(int i, double target);

extern float fp_alph_ini, fp_ffinc, fp_ffdec, fp_ffalph, fp_ffdtmax;
extern int fp_nfmin;

void relax_gloc()
{
    double dp = 0.0;
    if (atom.QC==0) {
      for (int i=1; i<=atom.natom; i++) {
	dp = dp + atom.vx[i]*atom.fx[i] + atom.vy[i]*atom.fy[i] + atom.vz[i]*atom.fz[i];
      }
      if (dp < 0.0) { vscale(0.0); }
    } else {
      for (int i=1; i<=atom.natom; i++) {
	if (atom.repatom[i]) {
	  dp = dp + atom.vx[i]*atom.fx[i] + atom.vy[i]*atom.fy[i] + atom.vz[i]*atom.fz[i];
	}
	// if (dp < 0.0) { vscale(i, 0.0); }
      }
      if (dp < 0.0) { vscale(0.0); }
    }
}

void relax_fire_reset()
{
  fire.sp = 0.0; fire.fnorm = 0.0; fire.vnorm = 0.0;
  fire.fire_alph = fp_alph_ini;
  fire.fire_alph_ini = fp_alph_ini;
  fire.nfmin = fp_nfmin;
  fire.ffinc = fp_ffinc; fire.ffdec = fp_ffdec; fire.ffalph = fp_ffalph;
  fire.ifire = 0;
  fire.ffdtmax = fp_ffdtmax*1.0e-15;
}

void relax_fire()
{
  /*
    Atom relaxation algorithm by Fast Inertial Relaxation Engine (FIRE)
    E.Bitzek, P.Koshinen, F.Gahler, M.Moseler and P.Gumbsch
    Physical Review Letters, 97, 170201 (2006) 
  */
  // Parameters for FIRE
  // Initial values are written in "myclass.h"
  //if (istep == 0) { fire.ffdtmax=dt*10.0; } // ROUGHLY ESTIMATED VALUE 
  //if (fire.ffdtmax < 1.0e-20) { fire.ffdtmax=dt*10.0; }

  //-----F1 : P=F*v 
  fire.ifire++;
  fire.sp=0.0;
  fire.fnorm=0.0;
  fire.vnorm=0.0;
  for (int i=1; i<=atom.natom; i++) {
    fire.sp = fire.sp + atom.fx[i]*atom.vx[i]+atom.fy[i]*atom.vy[i]+atom.fz[i]*atom.vz[i];
    fire.fnorm += atom.fx[i]*atom.fx[i]+atom.fy[i]*atom.fy[i]+atom.fz[i]*atom.fz[i];
    fire.vnorm += atom.vx[i]*atom.vx[i]+atom.vy[i]*atom.vy[i]+atom.vz[i]*atom.vz[i];
  }
  fire.fnorm=sqrt(fire.fnorm);
  fire.vnorm=sqrt(fire.vnorm);

  //-----F2 : v = (1-a)*v + a*F*|v|
  for (int i=1; i<=atom.natom; i++) {
    atom.vx[i]=(1.0-fire.fire_alph)*atom.vx[i]
      +fire.fire_alph*atom.fx[i]/fire.fnorm*fire.vnorm;
    atom.vy[i]=(1.0-fire.fire_alph)*atom.vy[i]
      +fire.fire_alph*atom.fy[i]/fire.fnorm*fire.vnorm;
    atom.vz[i]=(1.0-fire.fire_alph)*atom.vz[i]
      +fire.fire_alph*atom.fz[i]/fire.fnorm*fire.vnorm;
  }

  //-----F3 : P > 0
  if ((fire.sp > 0.0)&&(fire.ifire >= fire.nfmin)) {
    //dt = min(dt*ffinc,ffdtmax);
    dt = dt*fire.ffinc; if (dt > fire.ffdtmax) { dt = fire.ffdtmax; }
    fire.fire_alph = fire.fire_alph * fire.ffalph;
    //-----F4 : P <= 0
  }  else if (fire.sp <= 0.0) {
    fire.ifire = 0;
    dt = dt*fire.ffdec;
    fire.fire_alph = fire.fire_alph_ini;
    vscale(0.0);
  }
  // To avoid burst...
  double xx = atom.Enkin()*2.0/3.0/(double)atom.natom/1.380662e-23;
  if (xx > 500) {
    fire.ifire = 0;
    dt = dt*fire.ffdec;
    fire.fire_alph = fire.fire_alph_ini;
    vscale(0.0);
  }
  //      print *,'FIRE:', dt, ifire, fire_alph, ffdtmax
}
