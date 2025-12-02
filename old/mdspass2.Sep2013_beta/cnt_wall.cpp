#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"
#include <GL/glui.h>
#include <complex>

using namespace std;

void bond_set();
void calc_product(double x1, double y1, double z1, double x2, double y2, double z2,
		  double &xx, double &yy, double &zz, double &area);
void cnt_wall_set();
void cnt_wall_write(const char* fname); void cnt_wall_read(const char* fname);
void cnt_wall_setnormalvectors(); void cnt_wall_centering();
void cnt_wall_discard();

void cnt_wall_discard()
{
   // Deallocate arrays
  if (atom.wall_v) {
    for (int i=0; i<=atom.nwall; i++) { if (atom.wall_v[i]) { delete[] atom.wall_v[i]; } }
    delete[] atom.wall_v; atom.wall_v = NULL; }
  if (atom.wall_nvec) {
    for (int i=0; i<=atom.nwall; i++) { if (atom.wall_nvec[i]) { delete[] atom.wall_nvec[i]; } }
    delete[] atom.wall_nvec; atom.wall_nvec = NULL; }
  if (atom.wall_vrev) { delete[] atom.wall_vrev; atom.wall_vrev = NULL; }
  if (atom.wall_area) { delete[] atom.wall_area; atom.wall_area = NULL; }
  atom.nwall = 0;
}

void cnt_wall_set()
{
  std::cout<<"cnt_wall_set"<<std::endl;
  // Deallocate arrays
  if (atom.wall_v) {
    for (int i=0; i<=atom.nwall; i++) { if (atom.wall_v[i]) { delete[] atom.wall_v[i]; } }
    delete[] atom.wall_v; atom.wall_v = NULL; }
  if (atom.wall_nvec) {
    for (int i=0; i<=atom.nwall; i++) { if (atom.wall_nvec[i]) { delete[] atom.wall_nvec[i]; } }
    delete[] atom.wall_nvec; atom.wall_nvec = NULL; }
  if (atom.wall_vrev) { delete[] atom.wall_vrev; atom.wall_vrev = NULL; }
  if (atom.wall_area) { delete[] atom.wall_area; atom.wall_area = NULL; }
  //std::cout<<"cnt_wall_set deallocate done"<<std::endl;
  // Count number of walls
  bond_set(); //std::cout<<"Bond set."<<std::endl;
  int icnt = 0;
  for (int i=1; i<=atom.natom; i++) {
    if (atom.nneighbor[i] != 3) { continue;
    } else {
      icnt++; icnt++; icnt++;
    } }
  atom.nwall = icnt; printf("Number of CNT walls = %d\n",atom.nwall);
  // Allocate arrays
  atom.wall_vrev = new int[atom.nwall+1];
  atom.wall_area = new double[atom.nwall+1];
  for (int i=1; i<=atom.nwall; i++) { atom.wall_vrev[i] = 0; atom.wall_area[i] = 0; }
  atom.wall_v = new int*[atom.nwall+1];
  for (int i=0; i<=atom.nwall; i++) { atom.wall_v[i] = new int[3];
    atom.wall_v[i][0]=0; atom.wall_v[i][1]=0; atom.wall_v[i][2]=0; }
  atom.wall_nvec = new double*[atom.nwall+1];
  for (int i=0; i<=atom.nwall; i++) { atom.wall_nvec[i] = new double[3];
    atom.wall_nvec[i][0]=0; atom.wall_nvec[i][1]=0; atom.wall_nvec[i][2]=0; }
  icnt = 0;
  for (int i=1; i<=atom.natom; i++) {
    if (atom.nneighbor[i] != 3) { continue;
    } else {
      icnt++; if (icnt > atom.nwall) {std::cout<<"too many walls"<<std::endl;exit(0);}
      atom.wall_v[icnt][0] = i;
      atom.wall_v[icnt][1] = atom.neighbor[i][0];
      atom.wall_v[icnt][2] = atom.neighbor[i][1];
      icnt++; if (icnt > atom.nwall) {std::cout<<"too many walls"<<std::endl;exit(0);}
      atom.wall_v[icnt][0] = i;
      atom.wall_v[icnt][1] = atom.neighbor[i][1];
      atom.wall_v[icnt][2] = atom.neighbor[i][2];
      icnt++; if (icnt > atom.nwall) {std::cout<<"too many walls"<<std::endl;exit(0);}
      atom.wall_v[icnt][0] = i;
      atom.wall_v[icnt][1] = atom.neighbor[i][2];
      atom.wall_v[icnt][2] = atom.neighbor[i][0];
    } }
  
  // set-up normal vectors
  cnt_wall_setnormalvectors();
  // make normal vectors point inward
  cnt_wall_centering();
  std::cout<<"CNT wall set done."<<std::endl;
}

void calc_product(double x1, double y1, double z1, double x2, double y2, double z2,
		  double &xx, double &yy, double &zz, double &area)
{
  xx = y1*z2 - z1*y2;
  yy = z1*x2 - x1*z2;
  zz = x1*y2 - y1*x2;
  double norm = sqrt(xx*xx+yy*yy+zz*zz);
  xx=xx/norm; yy=yy/norm; zz=zz/norm;
  area = norm;
}

void cnt_wall_setnormalvectors()
{
  for (int i=1; i<=atom.nwall; i++) {
    double x1 = atom.rx[atom.wall_v[i][1]] - atom.rx[atom.wall_v[i][0]];
    double y1 = atom.ry[atom.wall_v[i][1]] - atom.ry[atom.wall_v[i][0]];
    double z1 = atom.rz[atom.wall_v[i][1]] - atom.rz[atom.wall_v[i][0]];
    double x2 = atom.rx[atom.wall_v[i][2]] - atom.rx[atom.wall_v[i][0]];
    double y2 = atom.ry[atom.wall_v[i][2]] - atom.ry[atom.wall_v[i][0]];
    double z2 = atom.rz[atom.wall_v[i][2]] - atom.rz[atom.wall_v[i][0]];
    double cube3 = cell.hmat[2][2];
    if (z1>cube3/2) {z1=z1-cube3;} if (z1<-cube3/2) {z1=z1+cube3;}
    if (z2>cube3/2) {z2=z2-cube3;} if (z2<-cube3/2) {z2=z2+cube3;}
    double xx, yy, zz, area;
    calc_product(x1,y1,z1,x2,y2,z2,xx,yy,zz,area);
    if (atom.wall_vrev[i] < 0) { xx=-xx; yy=-yy; zz=-zz; }
    atom.wall_nvec[i][0] = xx;
    atom.wall_nvec[i][1] = yy;
    atom.wall_nvec[i][2] = zz;
    atom.wall_area[i] = area/2.0;
  }
  //std::cout<<"Normal vectors of CNT walls have been set."<<std::endl;
}

void cnt_wall_centering()
{
  double xcnt=0.0, ycnt=0.0;
  for (int i=1; i<=atom.natom; i++) { xcnt=xcnt+atom.rx[i]; ycnt=ycnt+atom.ry[i]; }
  xcnt=xcnt/atom.natom; ycnt=ycnt/atom.natom;
  for (int i=1; i<=atom.nwall; i++) {
    double x0 = atom.rx[atom.wall_v[i][0]]; double y0 = atom.ry[atom.wall_v[i][0]];
    double x1 = atom.rx[atom.wall_v[i][1]]; double y1 = atom.ry[atom.wall_v[i][1]];
    double x2 = atom.rx[atom.wall_v[i][2]]; double y2 = atom.ry[atom.wall_v[i][2]];
    double cube3 = cell.hmat[2][2];
    double xc = (x0+x1+x2)/3-xcnt; double yc = (y0+y1+y2)/3-ycnt;
    double xv = atom.wall_nvec[i][0]; double yv = atom.wall_nvec[i][1];
    double ang_c, ang_v;
    if (xc == 0.0) {
      if (yc > 0.0) { ang_c = 90;
      } else { ang_c = 270; }
    } else {
      ang_c = atan(yc/xc)/M_PI*180;
      if ((xc>0)&&(yc<0)) { ang_c=ang_c+360; }
      if ((xc<0)&&(yc<0)) { ang_c=ang_c+180; }
      if ((xc<0)&&(yc>0)) { ang_c=ang_c+180; }
      if ((yc==0.0)&&(xc>0)) { ang_c=0; }
      if ((yc==0.0)&&(xc<0)) { ang_c=180; }
    }
    if (xv == 0.0) {
      if (yv > 0.0) { ang_v = 90;
      } else { ang_v = 270; }
    } else {
      ang_v = atan(yv/xv)/M_PI*180;
      if ((xv>0)&&(yv<0)) { ang_v=ang_v+360; }
      if ((xv<0)&&(yv<0)) { ang_v=ang_v+180; }
      if ((xv<0)&&(yv>0)) { ang_v=ang_v+180; }
      if ((yv==0.0)&&(xv>0)) { ang_v=0; }
      if ((yv==0.0)&&(xv<0)) { ang_v=180; }
    }
    atom.wall_vrev[i] = 1;
    if ((abs(ang_c-ang_v) < 20)||(abs(ang_c-ang_v) > 340)) {
      atom.wall_vrev[i] = -1;
      atom.wall_nvec[i][0]=-atom.wall_nvec[i][0];
      atom.wall_nvec[i][1]=-atom.wall_nvec[i][1];
      atom.wall_nvec[i][2]=-atom.wall_nvec[i][2]; }
  }
  //std::cout<<"And the vectors point to CNT center."<<std::endl;
}

void cnt_wall_read(const char* fname)
{
  printf ("Read CNT wall data from %s\n", fname); 
  std::ifstream fin(fname);
  // Deallocate arrays
  if (atom.wall_v) {
    for (int i=0; i<=atom.nwall; i++) { if (atom.wall_v[i]) { delete[] atom.wall_v[i]; } }
    delete[] atom.wall_v; atom.wall_v = NULL; }
  if (atom.wall_nvec) {
    for (int i=0; i<=atom.nwall; i++) { if (atom.wall_nvec[i]) { delete[] atom.wall_nvec[i]; } }
    delete[] atom.wall_nvec; atom.wall_nvec = NULL; }
  if (atom.wall_vrev) { delete[] atom.wall_vrev; atom.wall_vrev = NULL; }
  if (atom.wall_area) { delete[] atom.wall_area; atom.wall_area = NULL; }
  // Allocate arrays
  fin >> atom.nwall;
  printf ("Number of CNT walls = %d\n", atom.nwall); 
  atom.wall_vrev = new int[atom.nwall+1];
  atom.wall_area = new double[atom.nwall+1];
  atom.wall_v = new int*[atom.nwall+1];
  for (int i=0; i<=atom.nwall; i++) { atom.wall_v[i] = new int[3];
    atom.wall_v[i][0]=0; atom.wall_v[i][1]=0; atom.wall_v[i][2]=0; }
  atom.wall_nvec = new double*[atom.nwall+1];
  for (int i=0; i<=atom.nwall; i++) { atom.wall_nvec[i] = new double[3];
    atom.wall_nvec[i][0]=0; atom.wall_nvec[i][1]=0; atom.wall_nvec[i][2]=0; }
  //
  for (int i=1; i<=atom.nwall; i++) {
    fin >> atom.wall_v[i][0]; fin >> atom.wall_v[i][1]; fin >> atom.wall_v[i][2];
    fin >> atom.wall_vrev[i];
  }
  cnt_wall_setnormalvectors();
}

void cnt_wall_write(const char* fname)
{
   printf ("Write CNT wall data to %s\n", fname); 
   FILE *fp = fopen(fname, "w");
   fprintf(fp," %d\n",atom.nwall);
   for (int i=1; i<=atom.nwall; i++) {
     fprintf(fp," %d %d %d %d\n",atom.wall_v[i][0],atom.wall_v[i][1],atom.wall_v[i][2],
	     atom.wall_vrev[i]);
   }
   fclose(fp);
}
