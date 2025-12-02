#include <iostream>
#include<string.h>
#include<stdlib.h>
#include<fstream>
#include "myheader.h"

void adp_alloc();
void adp_param();
void mishin(double rr, double *param, double &val, double &grad);
void mishin_sc(double rr, double *param, double &val, double &grad);
void csw2(double rr, double *param, double &val, double &grad);
void csw2_sc(double rr, double *param, double &val, double &grad);
void poly_5(double rr, double *param, double &val, double &grad);
void exp_plus(double rr, double *param, double &val, double &grad);
void exp_plus_sc(double rr, double *param, double &val, double &grad);
double dsquare(double d);
void resetmat(double a[3][3]);
int atom_number(char* at);
int adptyp(char* at);

void e_force_adp()
{
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  //double rcut = 8.0e0; double rcut2 = rcut * rcut;
  int    self;
  int    h, k, l, typ1, typ2, uf, us, stresses, ipair;
  double value, grad, value_tail, grad_tail, grad_i, grad_j, p_sr_tail;
  double value_el, grad_el, ggrad_el;

  double phi_val, phi_grad, u_val, u_grad, w_val, w_grad, rho_val, rho_grad;
  double eb_val, eb_grad, rho_grad_j;

  // Initialization and setup
  if (adp.initialize) {
    adp_alloc();
    adp_param();
    adp.initialize = false;
  }
  rcut = adp.cut; rcut2 = rcut * rcut;

  // First loop (reset variables)
  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
    adp.rho[i] = 0.0; adp.gradF[i] = 0.0;
    adp.mu_x[i] = 0.0; adp.mu_y[i] = 0.0; adp.mu_z[i] = 0.0;
    adp.lambda_xx[i] = 0.0; adp.lambda_yy[i] = 0.0; adp.lambda_zz[i] = 0.0;
    adp.lambda_xy[i] = 0.0; adp.lambda_yz[i] = 0.0; adp.lambda_zx[i] = 0.0;
    for (int j=0; j<3; j++) { for (int k=0; k<3; k++) {
	atom.satom[i][j][k] = 0.0; } }
  }// FIRST LOOP END
  
  // Second loop (pair term and atomic density)
  for (int i=1; i<=atom.natom; i++) {
    typ1 = adptyp(atom.asp[i]);
    if (book.alistnum[i]>0) { for (int k=1; k<=book.alistnum[i]; k++) {
	typ2 = adptyp(atom.asp[j]);
	j  = book.alist[i][k][0];
	if (j<i) continue; //OK???
	ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	rr2 = atom.Dist2(i,j,ix,iy,iz)/ang/ang; //rr = sqrt(rr2);
	drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	drx /= ang; dry /= ang; drz /= ang;
	// pair potential
	if (rr2 < rcut2) { rr = sqrt(rr2);
	  self = 0; if (i==j) { self = 1; }
	  mishin_sc(rr, adp.param_pair, phi_val, phi_grad);
	  if (self) { phi_val *= 0.5; phi_grad *= 0.5; }
	  phi_grad /= rr;
	  atom.epot[i] += phi_val/2.0; atom.epot[j] += phi_val/2.0;
	  atom.fx[i] += drx*phi_grad; atom.fy[i] += dry*phi_grad; atom.fz[i] += drz*phi_grad;
	  atom.fx[j] -= drx*phi_grad; atom.fy[j] -= dry*phi_grad; atom.fz[j] -= drz*phi_grad;
	  atom.satom[i][0][0] += drx*drx*phi_grad/2.0;
	  atom.satom[i][0][1] += dry*drx*phi_grad/2.0;
	  atom.satom[i][0][2] += drz*drx*phi_grad/2.0;
	  atom.satom[i][1][1] += dry*dry*phi_grad/2.0;
	  atom.satom[i][1][2] += drz*dry*phi_grad/2.0;
	  atom.satom[i][2][2] += drz*drz*phi_grad/2.0;
	  atom.satom[j][0][0] += drx*drx*phi_grad/2.0;
	  atom.satom[j][0][1] += dry*drx*phi_grad/2.0;
	  atom.satom[j][0][2] += drz*drx*phi_grad/2.0;
	  atom.satom[j][1][1] += dry*dry*phi_grad/2.0;
	  atom.satom[j][1][2] += drz*dry*phi_grad/2.0;
	  atom.satom[j][2][2] += drz*drz*phi_grad/2.0;
	  // dipole distortion
	  exp_plus_sc(rr, adp.param_dp, u_val, u_grad);
	  if (self) { u_val *= 0.5; u_grad *= 0.5; }
	  //u_val /= rr; ###
	  adp.mu_x[i] += u_val*drx; adp.mu_x[j] -= u_val*drx;
	  adp.mu_y[i] += u_val*dry; adp.mu_y[j] -= u_val*dry;
	  adp.mu_z[i] += u_val*drz; adp.mu_z[j] -= u_val*drz;
	  // quadrupole distortion
	  exp_plus_sc(rr, adp.param_qp, w_val, w_grad);
	  if (self) { w_val *= 0.5; w_grad *= 0.5; }
	  //w_val /= rr2; ###
	  adp.lambda_xx[i] += w_val*drx*drx; adp.lambda_xx[j] += w_val*drx*drx;
	  adp.lambda_yy[i] += w_val*dry*dry; adp.lambda_yy[j] += w_val*dry*dry;
	  adp.lambda_zz[i] += w_val*drz*drz; adp.lambda_zz[j] += w_val*drz*drz;
	  adp.lambda_xy[i] += w_val*drx*dry; adp.lambda_xy[j] += w_val*drx*dry;
	  adp.lambda_yz[i] += w_val*dry*drz; adp.lambda_yz[j] += w_val*dry*drz;
	  adp.lambda_zx[i] += w_val*drz*drx; adp.lambda_zx[j] += w_val*drz*drx;
	  // atomic density
	  if (typ1 == typ2) {
	    csw2_sc(rr, adp.param_den, rho_val, rho_grad);
	    if (self) { rho_val *= 0.5; rho_grad *= 0.5; }
	    adp.rho[i] += rho_val; adp.rho[j] += rho_val;
	  } else {
	    csw2_sc(rr, adp.param_den, rho_val, rho_grad);
	    adp.rho[i] += rho_val; adp.rho[j] += rho_val;
	  }
	} //endif rr2<rcut2
      } } } // Second loop end
  for (int i=1; i<=atom.natom; i++) { // Second loop supplement
    // embedding energy
    poly_5(adp.rho[i], adp.param_eb, eb_val, eb_grad);
    atom.epot[i] += eb_val;
    adp.gradF[i] += eb_grad;
    // ADP energy
    double tmp = 0.0;
    tmp += dsquare(adp.mu_x[i]);
    tmp += dsquare(adp.mu_y[i]);
    tmp += dsquare(adp.mu_z[i]);
    adp.nu[i] = adp.lambda_xx[i] + adp.lambda_yy[i] + adp.lambda_zz[i];
    double tr = adp.nu[i] / 3.0;
    tmp += dsquare(adp.lambda_xx[i] - tr);
    tmp += dsquare(adp.lambda_yy[i] - tr);
    tmp += dsquare(adp.lambda_zz[i] - tr);
    tmp += dsquare(adp.lambda_xy[i])*2.0;
    tmp += dsquare(adp.lambda_yz[i])*2.0;
    tmp += dsquare(adp.lambda_zx[i])*2.0;
    tmp /= 2.0;
    atom.epot[i] += tmp;
  } // Second loop supplement end
  // Third loop (APD force)
  for (int i=1; i<=atom.natom; i++) {
    typ1 = adptyp(atom.asp[i]);
    if (book.alistnum[i]>0) { for (int k=1; k<=book.alistnum[i]; k++) {
	typ2 = adptyp(atom.asp[j]);
	j  = book.alist[i][k][0];
	if (j<i) continue; //OK???
	ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	rr2 = atom.Dist2(i,j,ix,iy,iz)/ang/ang; //rr = sqrt(rr2);
	drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	drx /= ang; dry /= ang; drz /= ang;
	if (rr2 < rcut2) { rr = sqrt(rr2);
	  self = 0; if (i==j) { self = 1; }
	  // rho
	  csw2_sc(rr, adp.param_den, rho_val, rho_grad);
	  if (typ1 == typ2) {
	    rho_grad_j = rho_grad;
	  } else {
	    csw2_sc(rr, adp.param_den, rho_val, rho_grad_j);
	  }
	  double tmp =
	    rho_grad * adp.gradF[i] + rho_grad_j * adp.gradF[j];
	  if (self) { tmp /= 2.0; }
	  tmp /= rr;
	  atom.fx[i] += drx*tmp; atom.fx[j] -= drx*tmp;
	  atom.fy[i] += dry*tmp; atom.fy[j] -= dry*tmp;
	  atom.fz[i] += drz*tmp; atom.fz[j] -= drz*tmp;
	  atom.satom[i][0][0] += drx*drx*tmp/2.0;
	  atom.satom[i][0][1] += dry*drx*tmp/2.0;
	  atom.satom[i][0][2] += drz*drx*tmp/2.0;
	  atom.satom[i][1][1] += dry*dry*tmp/2.0;
	  atom.satom[i][1][2] += drz*dry*tmp/2.0;
	  atom.satom[i][2][2] += drz*drz*tmp/2.0;
	  atom.satom[j][0][0] += drx*drx*tmp/2.0;
	  atom.satom[j][0][1] += dry*drx*tmp/2.0;
	  atom.satom[j][0][2] += drz*drx*tmp/2.0;
	  atom.satom[j][1][1] += dry*dry*tmp/2.0;
	  atom.satom[j][1][2] += drz*dry*tmp/2.0;
	  atom.satom[j][2][2] += drz*drz*tmp/2.0;
	  // dipole
	  double utmpx = adp.mu_x[i] - adp.mu_x[j];
	  double utmpy = adp.mu_y[i] - adp.mu_y[j];
	  double utmpz = adp.mu_z[i] - adp.mu_z[j];
	  //if (self) { utmpx /= 2.0; utmpy /= 2.0; utmpz /= 2.0; } // Bug in potfit? 
	  // (for small cell it does harm)
	  exp_plus_sc(rr, adp.param_dp, u_val, u_grad); // calculated again..
	  if (self) { u_val *= 0.5; u_grad *= 0.5; }
	  double utmp = (utmpx*drx + utmpy*dry + utmpz*drz) * u_grad;
	  utmp /= rr;
	  double tmpx = utmpx * u_val + utmp * drx;
	  double tmpy = utmpy * u_val + utmp * dry;
	  double tmpz = utmpz * u_val + utmp * drz;
	  atom.fx[i] += tmpx; atom.fx[j] -= tmpx;
	  atom.fy[i] += tmpy; atom.fy[j] -= tmpy;
	  atom.fz[i] += tmpz; atom.fz[j] -= tmpz;
	  atom.satom[i][0][0] += drx*tmpx/2.0;
	  atom.satom[i][0][1] += dry*tmpx/2.0;
	  atom.satom[i][0][2] += drz*tmpx/2.0;
	  atom.satom[i][1][1] += dry*tmpy/2.0;
	  atom.satom[i][1][2] += drz*tmpy/2.0;
	  atom.satom[i][2][2] += drz*tmpz/2.0;
	  atom.satom[j][0][0] += drx*tmpx/2.0;
	  atom.satom[j][0][1] += dry*tmpx/2.0;
	  atom.satom[j][0][2] += drz*tmpx/2.0;
	  atom.satom[j][1][1] += dry*tmpy/2.0;
	  atom.satom[j][1][2] += drz*tmpy/2.0;
	  atom.satom[j][2][2] += drz*tmpz/2.0;
	  // quadrupole
	  double wxx = adp.lambda_xx[i] + adp.lambda_xx[j];
	  double wyy = adp.lambda_yy[i] + adp.lambda_yy[j];
	  double wzz = adp.lambda_zz[i] + adp.lambda_zz[j];
	  double wxy = adp.lambda_xy[i] + adp.lambda_xy[j];
	  double wyz = adp.lambda_yz[i] + adp.lambda_yz[j];
	  double wzx = adp.lambda_zx[i] + adp.lambda_zx[j];
	  //if (self) { // Bug in potfit? (for small cell it does harm)
	  //  wxx /= 2.0; wyy /= 2.0; wzz /= 2.0;
	  //  wxy /= 2.0; wyz /= 2.0; wzx /= 2.0; }
	  double vx = (wxx*drx + wxy*dry + wzx*drz);
	  double vy = (wxy*drx + wyy*dry + wyz*drz);
	  double vz = (wzx*drx + wyz*dry + wzz*drz);
	  double nu = (adp.nu[i] + adp.nu[j])/3.0;
	  exp_plus_sc(rr, adp.param_qp, w_val, w_grad); // calculated again..
	  if (self) { w_val /= 2.0; w_grad /= 2.0; }
	  double f1 = w_val * 2.0;
	  double f2 = ((vx*drx+vy*dry+vz*drz)-nu*rr2)*w_grad-nu*f1*rr;
	  tmpx = f1 * vx + f2 * drx / rr;
	  tmpy = f1 * vy + f2 * dry / rr;
	  tmpz = f1 * vz + f2 * drz / rr;
	  atom.fx[i] += tmpx; atom.fx[j] -= tmpx;
	  atom.fy[i] += tmpy; atom.fy[j] -= tmpy;
	  atom.fz[i] += tmpz; atom.fz[j] -= tmpz;
	  atom.satom[i][0][0] += drx*tmpx/2.0;
	  atom.satom[i][0][1] += dry*tmpx/2.0;
	  atom.satom[i][0][2] += drz*tmpx/2.0;
	  atom.satom[i][1][1] += dry*tmpy/2.0;
	  atom.satom[i][1][2] += drz*tmpy/2.0;
	  atom.satom[i][2][2] += drz*tmpz/2.0;
	  atom.satom[j][0][0] += drx*tmpx/2.0;
	  atom.satom[j][0][1] += dry*tmpx/2.0;
	  atom.satom[j][0][2] += drz*tmpx/2.0;
	  atom.satom[j][1][1] += dry*tmpy/2.0;
	  atom.satom[j][1][2] += drz*tmpy/2.0;
	  atom.satom[j][2][2] += drz*tmpz/2.0;
	  

	} // rr2<rcut2
      } } } // Third loop end

 OUT:
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.epot[i] *= ev;
    atom.epotsum += atom.epot[i];
    atom.fx[i] *= ev/ang; atom.fy[i] *= ev/ang; atom.fz[i] *= ev/ang;
    //printf("%10d %20.12e %20.12e %20.12e\n",i-1,atom.fx[i]/ev*ang,atom.fy[i]/ev*ang,atom.fz[i]/ev*ang);
  }
    
  //Atomic stress
  for (int i=1; i<=atom.natom; i++) {
    atom.satom[i][1][0] = atom.satom[i][0][1];
    atom.satom[i][2][0] = atom.satom[i][0][2];
    atom.satom[i][2][1] = atom.satom[i][1][2]; }
  resetmat(cell.dmat);
  for (int i=1; i<=atom.natom; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	atom.satom[i][j][k] *= ev;
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

void mishin(double rr, double *param, double &val, double &grad)
{
  double z  = rr - param[3];
  double e  = exp(-param[5] * z);
  double p  = pow(z,param[4]);
  double pp = param[4] * p / z;
  double ep = -param[5] * e;
  val  = param[0] * p * e * (1.0 + param[1] * e) + param[2];
  grad = param[0] * ( pp  *  e * (1.0 + param[1] *  e)
		      + p * ep * (1.0 + param[1] *  e)
		      + p *  e * (      param[1] * ep) );
}
void mishin_sc(double rr, double *param, double &val, double &grad)
{
  double x = (rr - adp.cut)/param[6];
  double x2 = x*x; double x3 = x2*x; double x4 = x2*x2;
  double sc = x4 / (1.0 + x4); double scp = 4.0*x3/(1.0+x4)/(1.0+x4)/param[6];
  double z  = rr - param[3];
  //double e  = exp(-param[5] * z); //<-- should be this
  double e  = exp(-param[5] * rr); //<-- definition in potfit
  double p  = pow(z,param[4]);
  double pp = param[4] * p / z;
  double ep = -param[5] * e;
  double tmp  = param[0] * p * e * (1.0 + param[1] * e) + param[2];
  val  = tmp * sc;
  grad = param[0] * ( pp  *  e * (1.0 + param[1] *  e)
		      + p * ep * (1.0 + param[1] *  e)
		      + p *  e * (      param[1] * ep) ) * sc + tmp * scp;
}
void csw2(double rr, double *param, double &val, double &grad)
{
  double p  = pow(rr, param[3]);
  val  = ( 1.0 + param[0] * cos(param[1] * rr + param[2]) ) / p;
  grad = ( -param[0] * param[1] * sin(param[1] * rr + param[2]) ) / p
    - param[3] * val / rr;
}
void csw2_sc(double rr, double *param, double &val, double &grad)
{
  double x = (rr - adp.cut)/param[4];
  double x2 = x*x; double x3 = x2*x; double x4 = x2*x2;
  double sc = x4 / (1.0 + x4); double scp = 4.0*x3/(1.0+x4)/(1.0+x4)/param[4];
  double p  = pow(rr, param[3]);
  double val0 = ( 1.0 + param[0] * cos(param[1] * rr + param[2]) ) / p;
  //val  = ( 1.0 + param[0] * cos(param[1] * rr + param[2]) ) / p * sc;
  val = val0 * sc;
  grad = (( -param[0] * param[1] * sin(param[1] * rr + param[2]) ) / p
	  - param[3] * val0 / rr) * sc + val0 * scp;
}
void poly_5(double rr, double *param, double &val, double &grad)
{
  double  r1 = rr - 1.0;
  double dr1 = r1 * r1;
  val  = param[0] + 0.5 * param[1] * dr1
    + param[2]*r1*dr1 + param[3]*dr1*dr1 + param[4]*r1*dr1*dr1;
  grad = param[1] * r1
    + 3.0*param[2]*dr1 + 4.0*param[3]*r1*dr1 + 5.0*param[4]*dr1*dr1;
}
void exp_plus(double rr, double *param, double &val, double &grad)
{
  double tmp = param[0] * exp(-param[1]*rr);
  val  = tmp + param[2];
  grad = -param[1] * tmp;
}
void exp_plus_sc(double rr, double *param, double &val, double &grad)
{
  double x = (rr - adp.cut)/param[3];
  double x2 = x*x; double x3 = x2*x; double x4 = x2*x2;
  double sc = x4 / (1.0 + x4); double scp = 4.0*x3/(1.0+x4)/(1.0+x4)/param[3];
  double tmp = param[0] * exp(-param[1]*rr);
  val  = (tmp + param[2]) * sc;
  grad = (-param[1] * tmp) * sc + (tmp + param[2]) * scp;
}

void adp_param()
{
  adp.param_pair[0] = 0.00691025;
  adp.param_pair[1] = 29.92791121;
  adp.param_pair[2] = -0.19277842;
  adp.param_pair[3] = 0.00001165;
  adp.param_pair[4] = 0.98767802;
  adp.param_pair[5] = 0.15430821;
  adp.param_pair[6] = 0.5;
  adp.param_den[0] = -0.75198545;
  adp.param_den[1] = 0.69706211;
  adp.param_den[2] = -2.93829476;
  adp.param_den[3] = 1.36634379;
  adp.param_den[4] = 2.0;
  adp.param_eb[0] = -4.41231862;
  adp.param_eb[1] = 6.16520171;
  adp.param_eb[2] = -4.09634234;
  adp.param_eb[3] = -2.60663587;
  adp.param_eb[4] = 6.89527313;
  adp.param_dp[0] = 102.41749562;
  adp.param_dp[1] = 2.30668629;
  adp.param_dp[2] = -0.01064225;
  adp.param_dp[3] = 1.53956319;
  adp.param_qp[0] = -103.46596987;
  adp.param_qp[1] = 2.56364479;
  adp.param_qp[2] = 0.00128067;
  adp.param_qp[3] = 0.76080274;

}

void adp_alloc()
{
  if (adp.rho)        { delete[] adp.rho;       adp.rho       = NULL; }
  if (adp.mu_x)       { delete[] adp.mu_x;      adp.mu_x      = NULL; }
  if (adp.mu_y)       { delete[] adp.mu_y;      adp.mu_y      = NULL; }
  if (adp.mu_z)       { delete[] adp.mu_z;      adp.mu_z      = NULL; }
  if (adp.nu)         { delete[] adp.nu;        adp.nu        = NULL; }
  if (adp.gradF)      { delete[] adp.gradF;     adp.gradF     = NULL; }
  if (adp.lambda_xx)  { delete[] adp.lambda_xx; adp.lambda_xx = NULL; }
  if (adp.lambda_yy)  { delete[] adp.lambda_yy; adp.lambda_yy = NULL; }
  if (adp.lambda_zz)  { delete[] adp.lambda_zz; adp.lambda_zz = NULL; }
  if (adp.lambda_xy)  { delete[] adp.lambda_xy; adp.lambda_xy = NULL; }
  if (adp.lambda_yz)  { delete[] adp.lambda_yz; adp.lambda_yz = NULL; }
  if (adp.lambda_zx)  { delete[] adp.lambda_zx; adp.lambda_zx = NULL; }
  if (adp.param_pair) { delete[] adp.param_pair;adp.param_pair= NULL; }
  if (adp.param_den)  { delete[] adp.param_den; adp.param_den = NULL; }
  if (adp.param_eb)   { delete[] adp.param_eb;  adp.param_eb  = NULL; }
  if (adp.param_dp)   { delete[] adp.param_dp;  adp.param_dp  = NULL; }
  if (adp.param_qp)   { delete[] adp.param_qp;  adp.param_qp  = NULL; }
  adp.rho        = new double[atom.natom+1];
  adp.mu_x       = new double[atom.natom+1];
  adp.mu_y       = new double[atom.natom+1];
  adp.mu_z       = new double[atom.natom+1];
  adp.nu         = new double[atom.natom+1];
  adp.gradF      = new double[atom.natom+1];
  adp.lambda_xx  = new double[atom.natom+1];
  adp.lambda_yy  = new double[atom.natom+1];
  adp.lambda_zz  = new double[atom.natom+1];
  adp.lambda_xy  = new double[atom.natom+1];
  adp.lambda_yz  = new double[atom.natom+1];
  adp.lambda_zx  = new double[atom.natom+1];
  adp.param_pair = new double[7];
  adp.param_den  = new double[5];
  adp.param_eb   = new double[5];
  adp.param_dp   = new double[4];
  adp.param_qp   = new double[4];
}

int adptyp(char* at)
{
  return 1;

  int i = atom_number(at);
  if (i == 40) { return 1;
  } else if (i ==  8) { return 2;
  } else if (i == 39) { return 3;
  } else { return 0; }
}
