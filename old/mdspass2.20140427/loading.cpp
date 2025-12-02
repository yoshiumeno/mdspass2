#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

void cnt_wall_setnormalvectors(); void cnt_wall_centering();
extern int mode_cnt_corrugation, cnt_load_algo;
extern float cnt_ring_radius, cnt_ring_fmax, cnt_ring_sharpness;
extern float cnt_pressure, cnt_pressure_ftot, cnt_pressure_gpa;
extern float outermost_radius_f, cylinder_side_area_f;
extern float load_aux_ftot;
void loading();
void loading_old();
double repulsion(double r, double fmax, double sharpness);
void calc_center_xy(double &x, double &y);
void calc_cylinder_side_area(double center_x, double center_y, double &rr, double &area);
double theta(double x, double y);
void potential();

void loading()
{
  for (int i=1; i<=atom.natom; i++) { // Reset F by loading
    atom.fx_l[i] = 0; atom.fy_l[i] = 0; atom.fz_l[i] = 0;
  }

  if (mode_cnt_corrugation) { 
    double center_x, center_y, xx, yy, rr, ff, rr0, x0, y0;
    /*
    center_x = cell.hmat[0][0]/2.0;
    center_y = cell.hmat[1][1]/2.0;
    rr = 0;
    for (int i=1; i<=atom.natom; i++) {
      xx = atom.rx[i] - center_x;
      yy = atom.ry[i] - center_y;
      rr += sqrt(xx*xx+yy*yy);
    }
    rr /= (double)atom.natom;
    double cylinder_side_area = 2*M_PI*rr*cell.hmat[2][2]; //(in m^2)
    */
    double outermost_radius, cylinder_side_area;
    calc_center_xy(center_x, center_y);
    calc_cylinder_side_area(center_x, center_y, outermost_radius, cylinder_side_area);
    outermost_radius_f = (float)outermost_radius/ang;
    cylinder_side_area_f = (float)cylinder_side_area/ang/ang;
    if (cnt_load_algo == 0) { // Hydrostatic pressure (Wall) mode
      double factor0 = 1.0e-10*cnt_pressure;
      double area0 = 0.8731e-20;
      cnt_wall_setnormalvectors();
      cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.nwall; i++) {
	int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
	i0 = abs(i0); i1 = abs(i1); i2 = abs(i2);
	double factor = factor0*atom.wall_area[i]/area0;
	//factor = factor0; // <== If commented out, identical with old version
	double fxx = atom.wall_nvec[i][0]*factor;
	double fyy = atom.wall_nvec[i][1]*factor;
	double fzz = atom.wall_nvec[i][2]*factor;
	//cnt_pressure_ftot = cnt_pressure_ftot + sqrt(fxx*fxx+fyy*fyy+fzz*fzz); // Erroneous
	atom.fx_l[i0] += fxx; atom.fy_l[i0] += fyy; atom.fz_l[i0] += fzz;
	atom.fx_l[i1] += fxx; atom.fy_l[i1] += fyy; atom.fz_l[i1] += fzz;
	atom.fx_l[i2] += fxx; atom.fy_l[i2] += fyy; atom.fz_l[i2] += fzz;
      }
      //cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      //cnt_pressure_ftot *= 1e9; // total force (in nN)
      //printf("CNT pressure: Force per atom = %f (nN)\n",cnt_pressure_ftot);
    } else if (cnt_load_algo == 1) { // Ring-force mode
      double fmax, sharpness;
      cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.natom; i++) {
	center_x = cell.hmat[0][0]/2.0;
	center_y = cell.hmat[1][1]/2.0;
	xx = atom.rx[i] - center_x;
	yy = atom.ry[i] - center_y;
	rr = sqrt(xx*xx+yy*yy);
	rr0 = (double)cnt_ring_radius - rr*1.0e10;
	fmax = (double)cnt_ring_fmax;
	sharpness = (double)cnt_ring_sharpness;
	ff = repulsion(rr0, fmax, sharpness)*ev/ang;
	cnt_pressure_ftot = cnt_pressure_ftot + sqrt(ff*ff);
	atom.fx_l[i] += -ff * xx / rr;
	atom.fy_l[i] += -ff * yy / rr;
      }
      //cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      //cnt_pressure_ftot *= 1e9; // total force (in nN)
    } else if (cnt_load_algo == 2) { // Ring-rigid mode
      cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.natom; i++) {
	atom.fx_l[i] = atom.fx[i];
	atom.fy_l[i] = atom.fy[i];
	center_x = cell.hmat[0][0]/2.0;
	center_y = cell.hmat[1][1]/2.0;
	xx = atom.rx[i] - center_x;
	yy = atom.ry[i] - center_y;
	rr = sqrt(xx*xx+yy*yy);
	rr0 = (double)cnt_ring_radius - rr*1.0e10;
	if (rr0 < 0) {
	  x0 = atom.rx[i]; y0 = atom.ry[i];
	  atom.rx[i] = center_x + xx / rr * (double)cnt_ring_radius * ang;
	  atom.ry[i] = center_y + yy / rr * (double)cnt_ring_radius * ang;
	  //atom.vx[i] += (atom.rx[i] - x0) / dt;
	  //atom.vy[i] += (atom.ry[i] - y0) / dt;
	  atom.vx[i] = 0; atom.vy[i] = 0;
	}
      }
      potential();
      for (int i=1; i<=atom.natom; i++) {
	xx = atom.rx[i] - center_x;
	yy = atom.ry[i] - center_y;
	rr = sqrt(xx*xx+yy*yy);
	rr0 = (double)cnt_ring_radius - rr*1.0e10;
	if (rr0 <= 1.0e-5) {
	  x0 = atom.fx_l[i]; y0 = atom.fy_l[i];
	  atom.fx_l[i] = -atom.fx[i]; xx = atom.fx_l[i];
	  atom.fy_l[i] = -atom.fy[i]; yy = atom.fy_l[i];
	  //atom.fx[i] -= x0; atom.fy[i] -= y0;
	  cnt_pressure_ftot += sqrt(xx*xx+yy*yy);
	}
      }
      cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      cnt_pressure_ftot *= 1e9; // total force (in nN)
    }
    if (cnt_load_algo != 2) {
      // Bug fix 20140413
      cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.natom; i++) {
        double fxx = atom.fx_l[i]; double fyy = atom.fy_l[i]; double fzz = atom.fz_l[i]; 
        cnt_pressure_ftot += sqrt(fxx*fxx+fyy*fyy+fzz*fzz);
      }
      cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      cnt_pressure_ftot *= 1e9; // total force (in nN)
      //

      for (int i=1; i<=atom.natom; i++) { // Add F by loading
	atom.fx[i] += atom.fx_l[i]; atom.fy[i] += atom.fy_l[i]; atom.fz[i] += atom.fz_l[i];
      }
    }
  } // end of mode_cnt_corrugation
  
 // auxiliary constraint
  if (yzplane_punch) {
    yzplane_punch_ftot = 0.0;
    double center_x = cell.hmat[0][0]/2.0;
    double center_y = cell.hmat[1][1]/2.0;
    for (int i=1; i<=atom.natom; i++) {
      atom.fx_l[i] = atom.fx[i];
      double xx = atom.rx[i] - center_x;
      double xx0 = (double)yzplane_punch_d/2 - fabs(xx/ang);
      if (xx0 < 0) {
	double x0 = atom.rx[i];
	atom.rx[i] = center_x + (double)yzplane_punch_d/2 * ang * xx / fabs(xx);
	//atom.vx[i] += (atom.rx[i] - x0) / dt;
	atom.vx[i] = 0;
      }
    }
    potential();
    for (int i=1; i<=atom.natom; i++) {
      double xx = atom.rx[i] - center_x;
      double xx0 = (double)yzplane_punch_d/2 - fabs(xx/ang);
      if (xx0 <= 1.0e-5) {
	double x0 = atom.fx_l[i];
	atom.fx_l[i] = -atom.fx[i]; xx = atom.fx_l[i];
	yzplane_punch_ftot += sqrt(xx*xx);
      }
    }
    yzplane_punch_ftot *= 1e9; // total force (in nN)
  }    
  // end of auxiliary constraint
}

double repulsion(double r, double fmax, double sharpness)
{
  double f;
  //double fmax = 10.0;
  //double sharpness = 2.0;
  if (r < 0) {
    f = fmax; 
  } else {
    f = fmax*exp(-r*sharpness);
  }
  return f;
}

void loading_old()
{
 if (mode_cnt_corrugation) {
   double factor = 1.0e-10*cnt_pressure;
   cnt_wall_setnormalvectors();
   cnt_pressure_ftot = 0.0;
   for (int i=1; i<=atom.nwall; i++) {
     int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
     double fxx = atom.wall_nvec[i][0]*factor;
     double fyy = atom.wall_nvec[i][1]*factor;
     double fzz = atom.wall_nvec[i][2]*factor;
     cnt_pressure_ftot = cnt_pressure_ftot + fxx*fxx+fyy*fyy+fzz*fzz;
     atom.fx[i0] = atom.fx[i0] + fxx; atom.fy[i0] = atom.fy[i0] + fyy; atom.fz[i0] = atom.fz[i0] + fzz;
     atom.fx[i1] = atom.fx[i1] + fxx; atom.fy[i1] = atom.fy[i1] + fyy; atom.fz[i1] = atom.fz[i1] + fzz;
     atom.fx[i2] = atom.fx[i2] + fxx; atom.fy[i2] = atom.fy[i2] + fyy; atom.fz[i2] = atom.fz[i2] + fzz;
   }
   cnt_pressure_ftot = sqrt(cnt_pressure_ftot)/(double)atom.natom*1e9; // force per atom (in nN)
   //printf("CNT pressure: Force per atom = %f (nN)\n",cnt_pressure_ftot);
 }
}

void calc_center_xy(double &x, double &y)
{
  x = 0; y = 0;
  for (int i=1; i<=atom.natom; i++) {
    x += atom.rx[i]; y += atom.ry[i]; }
  x /= (double)atom.natom; y /= (double)atom.natom;
}

void calc_cylinder_side_area(double center_x, double center_y, double &rr, double &area)
{
  double d = 1.0*ang;
  double xx, yy, rr2, rr2_max;
  rr2_max = 0;
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i] - center_x;
    yy = atom.ry[i] - center_y;
    rr2 = xx*xx+yy*yy;
    if (rr2_max < rr2 ) { rr2_max = rr2; }
  }
  rr = 0;
  int cnt = 0;
  double rr_thres = sqrt(rr2_max) - d;
  double rr_thres2 = rr_thres * rr_thres;
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i] - center_x;
    yy = atom.ry[i] - center_y;
    rr2 = xx*xx+yy*yy;
    if (rr2 > rr_thres2) {
      rr += sqrt(rr2); cnt++;
    }
  }
  rr /= (double)cnt; // Calculated radius of outermost CNT (in m)
  //printf(" %d %f \n", cnt, rr/ang);
  area = 2.0*M_PI*rr*cell.hmat[2][2]; //(in m^2)
}
