#include <iostream>
#include "myheader.h"
double v(double rmeter);
double vp(double rmeter);

void e_force_morse()
{
  FLOAT rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  FLOAT hmat0[3], hmat1[3], hmat2[3];
  hmat0[0]=cell.hmat[0][0]; hmat1[0]=cell.hmat[1][0]; hmat2[0]=cell.hmat[2][0];
  hmat0[1]=cell.hmat[0][1]; hmat1[1]=cell.hmat[1][1]; hmat2[1]=cell.hmat[2][1];
  hmat0[2]=cell.hmat[0][2]; hmat1[2]=cell.hmat[1][2]; hmat2[2]=cell.hmat[2][2];

  for (int i=1; i<=atom.natom; i++) {
    atom.fx_float[i] = 0.0; atom.fy_float[i] = 0.0; atom.fz_float[i] = 0.0;  atom.epot_float[i] = 0.0;
    atom.rx_float[i] = atom.rx[i]; atom.ry_float[i] = atom.ry[i]; atom.rz_float[i] = atom.rz[i];
  }
  FLOAT ep = 0.3429; FLOAT al = 1.3588; FLOAT ro = 2.866; //FLOAT ev = 1.6021892e-19;
  //  FLOAT rcut = 8.0e-10;
  FLOAT rcut2 = rcut * rcut;

  for (int i=1; i<=atom.natom; i++) {
    if (book.alistnum[i]>0) {
      for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	//if (j>=i) {
	if ((j>=i)||(j<i)) { // <== Double counting algo
	  if ((atom.QC==0)||(atom.repatom[i]==1)||(atom.repatom[j]==1)) {
	    ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	    //rr2 = atom.Dist2(i,j,ix,iy,iz);
	    FLOAT xx = atom.rx_float[j] + hmat0[0]*ix+hmat0[1]*iy+hmat0[2]*iz - atom.rx_float[i];
	    FLOAT yy = atom.ry_float[j] + hmat1[0]*ix+hmat1[1]*iy+hmat1[2]*iz - atom.ry_float[i];
	    FLOAT zz = atom.rz_float[j] + hmat2[0]*ix+hmat2[1]*iy+hmat2[2]*iz - atom.rz_float[i];
	    rr2 = ( xx*xx + yy*yy + zz*zz );
	    if (rr2 < rcut2) {
	      rr = sqrt(rr2);
	      //drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	      drx = atom.rx_float[j] + hmat0[0]*ix+hmat0[1]*iy+hmat0[2]*iz - atom.rx_float[i];
	      dry = atom.ry_float[j] + hmat1[0]*ix+hmat1[1]*iy+hmat1[2]*iz - atom.ry_float[i];
	      drz = atom.rz_float[j] + hmat2[0]*ix+hmat2[1]*iy+hmat2[2]*iz - atom.rz_float[i];
	      //vp0 = vp(rr); v0 = v(rr);
	      double rang = rr * 1.0e10;
	      vp0 = -2.0*al*ep*(exp(-2.0*al*(rang-ro))-exp(-al*(rang-ro))) *ev*1.0e10;
	      v0 = ep*(exp(-2.0*al*(rang-ro))-2.0*exp(-al*(rang-ro))) *ev;

	      if ((atom.QC==0)||(atom.repatom[i]==1)) {
		atom.fx_float[i] = atom.fx_float[i]+vp0/rr*drx;
		atom.fy_float[i] = atom.fy_float[i]+vp0/rr*dry;
		atom.fz_float[i] = atom.fz_float[i]+vp0/rr*drz;
		atom.epot_float[i]=atom.epot_float[i]+v0/2.0;
		//if (i==j) { atom.epot_float[i]=atom.epot_float[i]+v0/2.0; } // <== Double counting algo
	      }
	    }
	  }
	} //j>=i
      }
    }
    //    printf("FAST %d %20.15e %20.15e\n",i,atom.fx[i],rr2);
  }
  //  if (istep==0) {exit(0);}

  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = atom.fx_float[i];
    atom.fy[i] = atom.fy_float[i];
    atom.fz[i] = atom.fz_float[i];
    atom.epot[i]=atom.epot_float[i];
  }

}


double v(double rmeter)
{
  double vval, rang, ep, al, ro;
  rang = rmeter * 1.0e10;
  //  ep = 0.2703;  al = 1.1646;  ro = 3.253;
  ep = 0.3429;  al = 1.3588;  ro = 2.866;
  //  ev = 1.6021892e-19;
  vval=ep*(exp(-2.0*al*(rang-ro))-2.0*exp(-al*(rang-ro))) *ev;
  return vval;
}

double vp(double rmeter)
{
  double vpval, rang, ep, al, ro;
  rang = rmeter * 1.0e10;
  //  ep = 0.2703;  al = 1.1646;  ro = 3.253;
  ep = 0.3429;  al = 1.3588;  ro = 2.866;
  //  ev = 1.6021892e-19;
  vpval=-2.0*al*ep*(exp(-2.0*al*(rang-ro))-exp(-al*(rang-ro))) *ev*1.0e10;
  return vpval;
}

