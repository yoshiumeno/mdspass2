#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include "myheader.h"
void tersoff_alloc();
void tersoff_setparam(int mparam);
double fcut(double r, double rr, double ss);
double fcutd(double r, double rr, double ss);
double v(double rmeter);
double vp(double rmeter);
void resetmat(double a[3][3]);

void e_force_tersoff()
{
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  rcut = 4.0e0; double rcut2 = rcut * rcut;

  double pi = M_PI, eps = 1.0e-20, ev2j = 1.6021892e-19, j2ev = 1.0/ev2j;
  int mode = 1; //1:38-14(B,C),2:39-8(Multi),1:Date's type
  int mparam = 0; //!0:Si(B), 1:Si(C) [see PRB 38 p.9902]
  //4:Date's potential(B),5:Yasukawa,6:Date's Potential(read from file)

  //if (strcmp(atom.potential_arg,"B")==0) mparam=0;
  //if (strcmp(atom.potential_arg,"C")==0) mparam=1;
  //if (strcmp(atom.potential_arg,"B2")==0) mparam=2;
  //if (strcmp(atom.potential_arg,"Sn")==0) mparam=3;

  if (tersoff.initialize) {
    std::cout<<"Initialize of the Tersoff potential..";
    tersoff_alloc();
    //printf("bookkeep algo = %d\n",book.algo);
    //if (book.algo != 2) { printf("Error in Tersoff: Bookkeep algo is not 2 (double)\n"); goto OUT; }
    if (atom.natom > 3000) { printf("ERROR: Please use TersoffNM instead.\n"); goto OUT; }
    tersoff_setparam(mparam);
    tersoff.initialize = false;
    std::cout<<"done"<<std::endl;
  }
 OUT:
  int istat = 1;
  double rc = tersoff.terss;
  rcut = rc * 1e10;
  double rc2 = rc*rc;

  int kns, kn, kns0;
  double rix, riy, riz, rjx, rjy, rjz, fij, dfij, rkx, rky, rkz;
  double dxij, dyij, dzij;

  //   Virial term reset
  cell.virx=0.0; cell.viry=0.0; cell.virz=0.0;
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
      cell.dmat[i][j] = 0.0;
      for (int ii=1; ii<=atom.natom; ii++) { atom.satom[ii][i][j] = 0.0; }  } }
  //========================================================
  //   Loop for b and z
  for (int i=1; i<=atom.natom; i++) {
    for (int j=1; j<=atom.natom; j++) {
      tersoff.zmat[i][j] = 0; tersoff.b[i][j] = 0; } }
  //kns = 0;
  for (int i=1; i<=atom.natom; i++) { // bloop
    rix = atom.rx[i]; riy = atom.ry[i]; riz = atom.rz[i];
    if (book.alistnum[i]>0) {
      for (int k0=1; k0<=book.alistnum[i]; k0++) { // 1200
	int j = book.alist[i][k0][0];
	int ix = book.alist[i][k0][1]; int iy = book.alist[i][k0][2]; int iz = book.alist[i][k0][3];
	double rij2 = atom.Dist2(i,j,ix,iy,iz);
	if ((rij2 < rc2)&&(rij2 > eps)) {
	  double rij = sqrt(rij2);
	  rjx = atom.rx[i] + atom.Dx(i,j,ix,iy,iz);
	  rjy = atom.ry[i] + atom.Dy(i,j,ix,iy,iz);
	  rjz = atom.rz[i] + atom.Dz(i,j,ix,iy,iz);
	  fij = fcut(rij,tersoff.terrr,tersoff.terss);
	  if (fij < eps) { continue; }
	  //--------------calc bij---------------------
	  for (int k1=1; k1<=book.alistnum[i]; k1++) {
	    int k = book.alist[i][k1][0];
	    ix = book.alist[i][k1][1]; iy = book.alist[i][k1][2]; iz = book.alist[i][k1][3];
	    rkx = atom.rx[i] + atom.Dx(i,k,ix,iy,iz);
	    rky = atom.ry[i] + atom.Dy(i,k,ix,iy,iz);
	    rkz = atom.rz[i] + atom.Dz(i,k,ix,iy,iz);
	    double rik2 = atom.Dist2(i,k,ix,iy,iz);
	    if (rik2 < eps) { continue; }
	    double rik = sqrt(rik2);
	    double fik = fcut(rik,tersoff.terrr,tersoff.terss);
	    if (fik < eps) { continue; }
	    double rjk2 = (rkx-rjx)*(rkx-rjx)+(rky-rjy)*(rky-rjy)+(rkz-rjz)*(rkz-rjz);
	    if (rjk2 < eps) { continue; }
	    double rjk = sqrt(rjk2);
	    double costh = (rik2+rij2-rjk2)/(2.0*rik*rij);
	    double g 
	      = 1.0+tersoff.terc2/tersoff.terd2
	      -tersoff.terc2/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh));
	    double zmatp;
	    if (mode == 1) {  zmatp = fik*g*exp(tersoff.termum*pow(rij-rik,tersoff.term));
	    } else { zmatp=fik*g*tersoff.terw; }
	    tersoff.zmat[i][j] = tersoff.zmat[i][j] + zmatp;
	  } // end of k1 loop
	  if (tersoff.b[i][j] != 0.0) {
	    printf("This Tersoff routine does not accept too small cells\n");
	    return; }
	  tersoff.b[i][j] 
	    = pow(1.0+pow(tersoff.terbeta*tersoff.zmat[i][j],tersoff.tern),-0.5/tersoff.tern);
	} // 1201
      } // 1200
    }
  } // bloop

  // Loop
  atom.epotsum = 0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
  }
  //kns = 0;
  for (int i=1; i<=atom.natom; i++) { // iloop
    rix = atom.rx[i]; riy = atom.ry[i]; riz = atom.rz[i];
    if (book.alistnum[i]>0) {
      for (int k0=1; k0<=book.alistnum[i]; k0++) { // 200
	int j = book.alist[i][k0][0];
	int ix = book.alist[i][k0][1]; int iy = book.alist[i][k0][2]; int iz = book.alist[i][k0][3];
	double rij2 = atom.Dist2(i,j,ix,iy,iz);
	if ((rij2 < rc2)&&(rij2 > eps)) {
	  double rij = sqrt(rij2);
	  rjx = atom.rx[i] + atom.Dx(i,j,ix,iy,iz);
	  rjy = atom.ry[i] + atom.Dy(i,j,ix,iy,iz);
	  rjz = atom.rz[i] + atom.Dz(i,j,ix,iy,iz);
	  fij = fcut(rij,tersoff.terrr,tersoff.terss);
	  dfij = fcutd(rij,tersoff.terrr,tersoff.terss);
	  if (fij < eps) { continue; }
	  double expmurij=exp(-tersoff.termu*rij);
	  double ushtij=-fij*tersoff.b[i][j]*tersoff.terbb*expmurij;
	  double urepij=fij*tersoff.teraa*exp(-tersoff.terlambda*rij);
	  double uij=urepij+ushtij;
	  atom.epot[i]=atom.epot[i]+uij*0.5;
	  // Atomic stress
	  double ad1 = 0.5*( (-tersoff.terlambda)*urepij + dfij/fij*urepij )/rij;
	  double adx=rjx-rix; double ady=rjy-riy; double adz=rjz-riz;
	  atom.satom[i][0][0] += ad1*adx*adx;
	  atom.satom[i][0][1] += ad1*ady*adx;
	  atom.satom[i][1][1] += ad1*ady*ady;
	  atom.satom[i][0][2] += ad1*adz*adx;
	  atom.satom[i][1][2] += ad1*adz*ady;
	  atom.satom[i][2][2] += ad1*adz*adz;
	  //
	  double dvdz=tersoff.terbb*fij*expmurij*tersoff.b[i][j]
	    *pow(tersoff.terbeta,tersoff.tern)*pow(tersoff.zmat[i][j],tersoff.tern)/2.0
	    /(1.0+pow(tersoff.terbeta*tersoff.zmat[i][j],tersoff.tern))/tersoff.zmat[i][j];
	  double ffx = 0.0; double ffy = 0.0; double ffz = 0.0;
	  for (int k1=1; k1<=book.alistnum[i]; k1++) { // kloop (300)
	    int k = book.alist[i][k1][0];
	    ix = book.alist[i][k1][1]; iy = book.alist[i][k1][2]; iz = book.alist[i][k1][3];
	    rkx = atom.rx[i] + atom.Dx(i,k,ix,iy,iz);
	    rky = atom.ry[i] + atom.Dy(i,k,ix,iy,iz);
	    rkz = atom.rz[i] + atom.Dz(i,k,ix,iy,iz);
	    double rik2 = atom.Dist2(i,k,ix,iy,iz);
	    if (rik2 < eps) { continue; }
	    double rik = sqrt(rik2);
	    double fik = fcut(rik,tersoff.terrr,tersoff.terss);
	    if (fik < eps) { continue; }
	    double rjk2 = (rkx-rjx)*(rkx-rjx)+(rky-rjy)*(rky-rjy)+(rkz-rjz)*(rkz-rjz);
	    if (rjk2 < eps) { continue; }
	    double rjk = sqrt(rjk2);
	    double costh = (rik2+rij2-rjk2)/(2.0*rik*rij);
	    double g 
	      = 1.0+tersoff.terc2/tersoff.terd2
	      -tersoff.terc2/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh));
	    double dzdr1, dzdc;
	    if (mode == 1) {
	      dzdr1=(double)tersoff.term*tersoff.termum
		*pow(rij-rik,tersoff.term-1)
		*fik*exp(tersoff.termum*pow(rij-rik,tersoff.term))
		*(1.0+tersoff.terc2/tersoff.terd2
		  -tersoff.terc2/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh)));
	      dzdc=fik*exp(tersoff.termum*pow(rij-rik,tersoff.term))
		*(-2.0*tersoff.terc2*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh));
	    } else {
	      dzdr1=0.0;
	      dzdc=fik*tersoff.terw
		*(-2.0*tersoff.terc2*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh)); }
	    double dzdr2=fcutd(rik,tersoff.terrr,tersoff.terss)
	      *exp(tersoff.termum*pow(rij-rik,tersoff.term))
	      *(1.0+tersoff.terc2/tersoff.terd2
		-tersoff.terc2/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh)))
	      -dzdr1;
	    dxij=rix-rjx; dyij=riy-rjy; dzij=riz-rjz;
	    double dxki=rkx-rix; double dyki=rky-riy; double dzki=rkz-riz;
	    double dcdrx=1.0/rik*(dxij/rij+dxki*costh/rik)
	      -1.0/rij*(dxki/rik+dxij*costh/rij);
	    double dcdry=1.0/rik*(dyij/rij+dyki*costh/rik)
	      -1.0/rij*(dyki/rik+dyij*costh/rij);
	    double dcdrz=1.0/rik*(dzij/rij+dzki*costh/rik)
	      -1.0/rij*(dzki/rik+dzij*costh/rij);
	    ffx=ffx-dvdz*(dzdr1*dxij/rij-dzdr2*dxki/rik+dzdc*dcdrx)/2.0;
	    ffy=ffy-dvdz*(dzdr1*dyij/rij-dzdr2*dyki/rik+dzdc*dcdry)/2.0;
	    ffz=ffz-dvdz*(dzdr1*dzij/rij-dzdr2*dzki/rik+dzdc*dcdrz)/2.0;
	    // Atomic stress
	    double ad1 = 0.5*dvdz*dzdr1/rij;
	    double ad2 = 0.5*dvdz*dzdr2/rik;
	    atom.satom[i][0][0] += ad1*dxij*dxij+ad2*dxki*dxki;
	    atom.satom[i][0][1] += ad1*dyij*dxij+ad2*dyki*dxki;
	    atom.satom[i][1][1] += ad1*dyij*dyij+ad2*dyki*dyki;
	    atom.satom[i][0][2] += ad1*dzij*dxij+ad2*dzki*dxki;
	    atom.satom[i][1][2] += ad1*dzij*dyij+ad2*dzki*dyki;
	    atom.satom[i][2][2] += ad1*dzij*dzij+ad2*dzki*dzki;
	    double dxjk=rkx-rjx; double dyjk=rky-rjy; double dzjk=rkz-rjz;
	    double ad0 = 0.5*dvdz*dzdc;
	    ad1 = -1.0/rik/rij;
	    ad2 = 1.0/rik/rij-costh/rik2;
	    double ad3 = 1.0/rik/rij-costh/rij2;
	    atom.satom[i][0][0] += ad0*(ad1*dxjk*dxjk+ad2*dxki*dxki+ad3*dxij*dxij);
	    atom.satom[i][0][1] += ad0*(ad1*dyjk*dxjk+ad2*dyki*dxki+ad3*dyij*dxij);
	    atom.satom[i][1][1] += ad0*(ad1*dyjk*dyjk+ad2*dyki*dyki+ad3*dyij*dyij);
	    atom.satom[i][0][2] += ad0*(ad1*dzjk*dxjk+ad2*dzki*dxki+ad3*dzij*dxij);
	    atom.satom[i][1][2] += ad0*(ad1*dzjk*dyjk+ad2*dzki*dyki+ad3*dzij*dyij);
	    atom.satom[i][2][2] += ad0*(ad1*dzjk*dzjk+ad2*dzki*dzki+ad3*dzij*dzij);
	    //
	  } // kloop (300)
	  double dvdr1=tersoff.teraa*exp(-tersoff.terlambda*rij)
	    *(dfij-tersoff.terlambda*fij)-tersoff.b[i][j]*tersoff.terbb*expmurij
	    *(dfij-tersoff.termu*fij);
	  double dvdr2=tersoff.teraa*exp(-tersoff.terlambda*rij)
	    *(dfij-tersoff.terlambda*fij)-tersoff.b[j][i]*tersoff.terbb*expmurij
	    *(dfij-tersoff.termu*fij);
	  atom.fx[i]=atom.fx[i]-((dvdr1+dvdr2)*dxij/rij)/2.0;
	  atom.fy[i]=atom.fy[i]-((dvdr1+dvdr2)*dyij/rij)/2.0;
	  atom.fz[i]=atom.fz[i]-((dvdr1+dvdr2)*dzij/rij)/2.0;
	  
	  atom.fx[i]=atom.fx[i]+ffx;
	  atom.fy[i]=atom.fy[i]+ffy;
	  atom.fz[i]=atom.fz[i]+ffz;
	  // Atomic stress
	  ad1 = 0.5*(fij*tersoff.termu-dfij)*tersoff.terbb*
	    (tersoff.b[i][j]+tersoff.b[j][i])/2.0*expmurij/rij;
	  adx=rjx-rix; ady=rjy-riy; adz=rjz-riz;
	  atom.satom[i][0][0] += ad1*adx*adx;
	  atom.satom[i][0][1] += ad1*ady*adx;
	  atom.satom[i][1][1] += ad1*ady*ady;
	  atom.satom[i][0][2] += ad1*adz*adx;
	  atom.satom[i][1][2] += ad1*adz*ady;
	  atom.satom[i][2][2] += ad1*adz*adz;
	  //
	} // 201
      } // 200
    }
  } // iloop
  
  // u_ji/r_i
  //kns = 0;
  for (int i=1; i<=atom.natom; i++) { // iloop2
    rix = atom.rx[i]; riy = atom.ry[i]; riz = atom.rz[i];
    if (book.alistnum[i]>0) {
      for (int k0=1; k0<=book.alistnum[i]; k0++) { // 2200
	int j = book.alist[i][k0][0];
	int ix = book.alist[i][k0][1]; int iy = book.alist[i][k0][2]; int iz = book.alist[i][k0][3];  
	double rij2 = atom.Dist2(i,j,ix,iy,iz);
	if ((rij2 < rc2)&&(rij2 > eps)) {
	  double rij = sqrt(rij2);
	  rjx = atom.rx[i] + atom.Dx(i,j,ix,iy,iz);
	  rjy = atom.ry[i] + atom.Dy(i,j,ix,iy,iz);
	  rjz = atom.rz[i] + atom.Dz(i,j,ix,iy,iz);
	  fij = fcut(rij,tersoff.terrr,tersoff.terss);
	  dfij = fcutd(rij,tersoff.terrr,tersoff.terss);
	  if (fij < eps) { continue; }
	  double ffx = 0.0; double ffy = 0.0; double ffz = 0.0;
	  double dvdz0=tersoff.terbb*fij*exp(-tersoff.termu*rij)*tersoff.b[j][i]
	    *pow(tersoff.terbeta,tersoff.tern)*pow(tersoff.zmat[j][i],tersoff.tern)/2.0
	    /(1.0+pow(tersoff.terbeta*tersoff.zmat[j][i],tersoff.tern))/tersoff.zmat[j][i];
	  for (int k1=1; k1<=book.alistnum[j]; k1++) { // kloop (300)
	    int k = book.alist[j][k1][0];
	    ix = book.alist[j][k1][1]; iy = book.alist[j][k1][2]; iz = book.alist[j][k1][3];
	    rkx = rjx + atom.Dx(j,k,ix,iy,iz);
	    rky = rjy + atom.Dy(j,k,ix,iy,iz);
	    rkz = rjz + atom.Dz(j,k,ix,iy,iz);
	    double rik2 = (rkx-rix)*(rkx-rix)+(rky-riy)*(rky-riy)+(rkz-riz)*(rkz-riz);
	    if (rik2 < eps) { continue; }
	    double rik = sqrt(rik2);
	    double rjk2 = (rkx-rjx)*(rkx-rjx)+(rky-rjy)*(rky-rjy)+(rkz-rjz)*(rkz-rjz);
	    if (rjk2 < eps) { continue; }
	    double rjk = sqrt(rjk2);
	    double fjk = fcut(rjk,tersoff.terrr,tersoff.terss);
	    if (fjk < eps) { continue; }
	    // u_ij/r_i
	    double costh=(rjk2+rij2-rik2)/(2.0*rjk*rij);
	    double g
	      = 1.0+tersoff.terc2/tersoff.terd2
	      -tersoff.terc2/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh));
	    double dzdr, dzdc;
	    if (mode == 1) {
	      dzdr=tersoff.termum*(double)tersoff.term
		*pow(rij-rjk,tersoff.term-1)
		*fjk*exp(tersoff.termum*pow(rij-rjk,tersoff.term))
	 	*(1.0+tersoff.terc2/tersoff.terd2
		  -tersoff.terc2/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh)));
	      dzdc=fjk*exp(tersoff.termum*pow(rij-rjk,tersoff.term))
		*(-2.0*tersoff.terc2*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh));
	    } else {
	      dzdr=0.0;
	      dzdc=fjk*tersoff.terw
		*(-2.0*tersoff.terc2*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh)); }
	    double ad1=(rjk2-rik2-rij2)/(2.0*rij2*rij*rjk);
	    double dcdrx=ad1*(rjx-rix)+1.0/(rij*rjk)*(rkx-rix);
	    double dcdry=ad1*(rjy-riy)+1.0/(rij*rjk)*(rky-riy);
	    double dcdrz=ad1*(rjz-riz)+1.0/(rij*rjk)*(rkz-riz);
	    dxij=rix-rjx;dyij=riy-rjy;dzij=riz-rjz;
	    ffx=ffx-dvdz0*(dzdr*dxij/rij+dzdc*dcdrx)/2.0;
	    ffy=ffy-dvdz0*(dzdr*dyij/rij+dzdc*dcdry)/2.0;
	    ffz=ffz-dvdz0*(dzdr*dzij/rij+dzdc*dcdrz)/2.0;
	    //u_jk/r_i
	    double dvdz=tersoff.terbb*fjk*exp(-tersoff.termu*rjk)*tersoff.b[j][k]
	      *pow(tersoff.terbeta,tersoff.tern)*pow(tersoff.zmat[j][k],tersoff.tern)/2.0
	      /(1.0+pow(tersoff.terbeta*tersoff.zmat[j][k],tersoff.tern))/tersoff.zmat[j][k];
	    if (mode == 1) {
	      dzdr=-(double)tersoff.term*tersoff.termum
		*pow(rjk-rij,tersoff.term-1)
		*fij*g*exp(tersoff.termum*pow(rjk-rij,tersoff.term))
		+dfij*exp(tersoff.termum*pow(rjk-rij,tersoff.term))*g;
	      dzdc=fij*exp(tersoff.termum*pow(rjk-rij,tersoff.term))
		*(-2.0*tersoff.terc2*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh));
	    } else {
	      dzdr=0.0;
	      dzdc=fij*tersoff.terw
		*(-2.0*tersoff.terc2*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh))
		/(tersoff.terd2+(tersoff.terh-costh)*(tersoff.terh-costh)); }
	    double ad5=(rjk2-rik2-rij2)/2.0/rij2/rij/rjk;
	    double ad6=1.0/rij/rjk;
	    dcdrx=(ad5*(rjx-rix)+ad6*(rkx-rix));
	    dcdry=(ad5*(rjy-riy)+ad6*(rky-riy));
	    dcdrz=(ad5*(rjz-riz)+ad6*(rkz-riz));
	    dxij=rix-rjx;dyij=riy-rjy;dzij=riz-rjz;
	    ffx=ffx-dvdz*(dzdr*dxij/rij+dzdc*dcdrx)/2.0;
	    ffy=ffy-dvdz*(dzdr*dyij/rij+dzdc*dcdry)/2.0;
	    ffz=ffz-dvdz*(dzdr*dzij/rij+dzdc*dcdrz)/2.0;
	  } // 2300
	  atom.fx[i]=atom.fx[i]+ffx;
	  atom.fy[i]=atom.fy[i]+ffy;
	  atom.fz[i]=atom.fz[i]+ffz;
	} // 2201
      } //2200
    }
  } // iloop2

  /*
  for (int i=1; i<=atom.natom; i++) {
    if ((atom.QC==0)||(atom.repatom[i]==1)) {
      if (book.alistnum[i]>0) {
	for (int k=1; k<=book.alistnum[i]; k++) {
	  j  = book.alist[i][k][0];
	  ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	  rr2 = atom.Dist2(i,j,ix,iy,iz);
	  if (rr2 < rcut2) {
	    rr = sqrt(rr2);
	    drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	    atom.fx[i] = atom.fx[i]+vp(rr)/rr*drx;
	    atom.fy[i] = atom.fy[i]+vp(rr)/rr*dry;
	    atom.fz[i] = atom.fz[i]+vp(rr)/rr*drz;
	    atom.epot[i]=atom.epot[i]+v(rr)/2.0;
	  }
	}
      }
    }
  }
  */
  //Etot
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.epotsum += atom.epot[i]; }
  //Atomic stress
  for (int i=1; i<=atom.natom; i++) {
    atom.satom[i][1][0] = atom.satom[i][0][1];
    atom.satom[i][2][0] = atom.satom[i][0][2];
    atom.satom[i][2][1] = atom.satom[i][1][2]; }
  resetmat(cell.dmat);
  for (int i=1; i<=atom.natom; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	cell.dmat[j][k] -= atom.satom[i][j][k]; } } }
  cell.virx = cell.dmat[0][0];
  cell.viry = cell.dmat[1][1];
  cell.virz = cell.dmat[2][2];
  cell.volume = cell.Getvolume();
  for (int i=1; i<=atom.natom; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	atom.satom[i][j][k] *= (double)atom.natom / cell.volume; } } }
}


void tersoff_alloc()
{
  tersoff.zmat = new double*[atom.natom+1];
  tersoff.b = new double*[atom.natom+1];
  for (int i=0; i<atom.natom+1; i++) {
    tersoff.zmat[i] = new double[atom.natom+1];
    tersoff.b[i] = new double[atom.natom+1]; }
}

void tersoff_setparam(int mparam)
{
  double pi = M_PI, eps = 1.0e-20, ev2j = 1.6021892e-19, j2ev = 1.0/ev2j;
  if (mparam == 0) {     // Si(B) (for surface structure)
    tersoff.terrr=2.80e-10;
    tersoff.terss=3.20e-10;
    tersoff.teraa=3.2647e3*ev2j;
    tersoff.terbb=9.5373e1*ev2j;
    tersoff.terlambda=3.2394e10;
    tersoff.termu=1.3258e10;
    tersoff.terbeta=3.3675e-1;
    tersoff.tern=2.2956e1;
    tersoff.term=3;
    tersoff.terc=4.8381e0;
    tersoff.terd=2.0417e0;
    tersoff.terh=0.0e0;
    tersoff.terw=1.0e0;
  } else if (mparam == 1) {  // Si(C) (for elastic property)
    tersoff.terrr=2.70e-10;
    tersoff.terss=3.00e-10;
    tersoff.teraa=1.8308e3*ev2j;
    tersoff.terbb=4.7118e2*ev2j;
    tersoff.terlambda=2.4799e10;
    tersoff.termu=1.7322e10;
    tersoff.terbeta=1.0999e-6;
    tersoff.tern=7.8734e-1;
    tersoff.term=3;
    tersoff.terc=1.0039e5;
    tersoff.terd=1.6218e1;
    tersoff.terh=-5.9826e-1;
    tersoff.terw=1.0e0;
  } else if (mparam == 2) { // Si(B*) (for phase transformation)
    tersoff.terrr=2.65e-10;
    tersoff.terss=2.85e-10;
    tersoff.teraa=3.2647e3*ev2j;
    tersoff.terbb=9.5373e1*ev2j;
    tersoff.terlambda=3.2394e10;
    tersoff.termu=1.3258e10;
    tersoff.terbeta=3.3675e-1;
    tersoff.tern=2.2956e1;
    tersoff.term=3;
    tersoff.terc=4.8381e0;
    tersoff.terd=2.0417e0;
    tersoff.terh=0.0e0;
    tersoff.terw=1.0e0;
  } else if (mparam == 3) { // Sn (Negami)
    tersoff.terrr=4.15e-10;
    tersoff.terss=4.45e-10;
    tersoff.teraa=1232.0e0*ev2j;
    tersoff.terbb=130.0e0*ev2j;
    tersoff.terlambda=2.277e10;
    tersoff.termu=1.203e10;
    tersoff.tern=0.9907e0;
    tersoff.terbeta=pow(0.1590e0,(1.0/tersoff.tern));
    tersoff.term=1;
    tersoff.terc=312.78e0;
    tersoff.terd=14.489e0;
    tersoff.terh=-0.4604e0;
    tersoff.terw=1.0e0;
  } else {
    printf("Please specify mparam\n");exit(0);
  }
  /*
!     FOR CARBON
!      terrr(2)=1.80e-10
!      terss(2)=2.10e-10
!      teraa(2)=1.3936e3*ev2j
!      terbb(2)=3.467e2*ev2j
!      terlambda(2)=3.4879e10
!      termu(2)=2.2119e10
!      terbeta(2)=1.5724e-7
!      tern(2)=7.2751e-1
!      term(2)=3
!      terc(2)=3.8049e4
!      terd(2)=4.384e0
!      terh(2)=-5.7058e-1
!      terw(2)=1.0e0
  */
  tersoff.terc2=tersoff.terc*tersoff.terc;
  tersoff.terd2=tersoff.terd*tersoff.terd;
  tersoff.termum=pow(tersoff.termu,tersoff.term);
}

double fcut(double r, double rr, double ss)
{
  double pi = M_PI, f;
  if (r < rr) {
    f=1.0;
  } else if (r <= ss) {
    f=0.5+0.5*cos(pi*(r-rr)/(ss-rr));
  } else {
    f=0.0;
  }
  return f;
}
double fcutd(double r, double rr, double ss)
{
  double pi = M_PI, fd;
  if ((r < ss)&&(r > rr)) {
    fd=-pi/2.0*sin(pi*(r-rr)/(ss-rr))/(ss-rr);
  } else {
    fd=0.0;
  }
  return fd;
}
