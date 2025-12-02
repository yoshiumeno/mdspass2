#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "myheader.h"
void matcpy(double a[3][3], double b[3][3]);
void resetmat(double a[3][3]);
void resetmat6(double a[6][6]);
void potential();
void e_force_dipole(int mode);

void force_check(double dis)
{
  double threshold=30.0;             // Error exceeding this value is shown
  double rxk, ryk, rzk, fx_cal, fy_cal, fz_cal, enpot0, fx_ene, fy_ene, fz_ene;
  double errx, erry, errz;
  dis *= ang;

  printf("Force check ...\n");
  double xx0,xx1;

  /*
  rxk=atom.rx[10];
  for (int i=-10; i<=10; i++) {
    atom.rx[10]=rxk+dis*(double)i;
    e_force_dipole(0);
    atom.rx[10]=rxk;
    e_force_dipole(1);
    printf("%d %20.10e\n",i,atom.epotsum/ev);
  }
  exit(0);
  */

  for (int i=1; i<=atom.natom; i++) {
    rxk=atom.rx[i]; ryk=atom.ry[i]; rzk=atom.rz[i];
    potential();
    //e_force_dipole(0);
    fx_cal=atom.fx[i]; fy_cal=atom.fy[i]; fz_cal=atom.fz[i];

    atom.rx[i]=rxk-dis;
    potential();
    //e_force_dipole(1);
    enpot0=atom.epotsum;
    atom.rx[i]=rxk+dis;
    potential();
    //e_force_dipole(1);
    //if (i==1) { xx1=dipole.val1; }
    //if (i==1) { printf("%e %e\n",-(xx1-xx0)/(dis*2/ang),dipole.val1); }
    //if (i==1) { printf("%e %e\n",xx1,xx0); }
    fx_ene=-(atom.epotsum-enpot0)/(dis*2.0);
    atom.rx[i]=rxk;
    errx=fabs(fx_cal/fx_ene*100.0-100);
    if (fx_ene==0.0) { errx=-1.0; }
    
    atom.ry[i]=ryk-dis;
    potential();
    //e_force_dipole(1);
    enpot0=atom.epotsum;
    atom.ry[i]=ryk+dis;
    potential();
    //e_force_dipole(1);
    fy_ene=-(atom.epotsum-enpot0)/(dis*2.0);
    atom.ry[i]=ryk;
    erry=fabs(fy_cal/fy_ene*100.0-100);
    if (fy_ene==0.0) { erry=-1.0; }

    atom.rz[i]=rzk-dis;
    potential();
    //e_force_dipole(1);
    enpot0=atom.epotsum;
    atom.rz[i]=rzk+dis;
    potential();
    //e_force_dipole(1);
    fz_ene=-(atom.epotsum-enpot0)/(dis*2.0);
    atom.rz[i]=rzk;
    errz=fabs(fz_cal/fz_ene*100.0-100);
    if (fz_ene==0.0) { errz=-1.0; }

    fx_cal /= ev/ang; fy_cal /= ev/ang; fz_cal /= ev/ang;
    fx_ene /= ev/ang; fy_ene /= ev/ang; fz_ene /= ev/ang;
    printf("Atom = %4d (x) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f %10.4f (%)\n",
	   i,fx_ene,fx_cal,fabs(fx_ene-fx_cal),errx);
    printf("            (y) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f %10.4f (%)\n",
	   fy_ene,fy_cal,fabs(fy_ene-fy_cal),erry);
    printf("            (z) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f %10.4f (%)\n",
	   fz_ene,fz_cal,fabs(fz_ene-fz_cal),errz);
  }
}

