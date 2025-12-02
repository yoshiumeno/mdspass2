#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"
#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

//double v(double rmeter);
//double vp(double rmeter);
//void e_force_morse_nobook();
//void e_force_morse();
void inverse(double mat[3][3], double imat[3][3]);
void potential();
void vscale();
void vscale(double target);
void vscale(int i, double target);
void writedata();
void stretch(double x, double y, double z);
void loading();
void velocity_random(double target);
void relax_gloc(); void relax_fire();
//void writeconfig(const char* fname);
//void capture();
void velocity_verlet_a();
void velocity_verlet_b();
void gear_pc_predic();
void gear_pc_correc();
void gear_pc_nphpr_predic();
void gear_pc_nphpr_correc();
void velocity_random(double target);
void recipe();
extern float ex,dexdt;
extern float ey,deydt;
extern float ez,dezdt;
extern int relax_algo;
extern float dtm;
//extern int confwrint, autocap;
extern int itolfor; extern float tolfor;
extern int itolstep; extern int tolstep;
extern int irecipe;
extern bool iprdamper, ivscale, irelax;
extern int istep0;
//void writedata(FILE *fp);

void md_set()
{
  integral_a[0] = velocity_verlet_a;
  integral_b[0] = velocity_verlet_b;
  integral_a[1] = gear_pc_nphpr_predic;
  integral_b[1] = gear_pc_nphpr_correc;
  //integral_b[1] = gear_pc_correc;
  if ((ensemble == 1)||(ensemble == 4)||(ensemble == 6)) { 
    velocity_random(temp_set); ivscale = true;
  } else {
    ivscale = false;
  }
  if (ensemble <= 2) {
    integral_type = 0;
  } else if (ensemble <= 7) {
    integral_type = 1;
  }
  if ((ensemble == 2)||(ensemble == 7)) {
    irelax = true;
  } else {
    irelax = false;
  }
  if ((ensemble >= 5)&&(ensemble <=7)) {
    iprdamper = true;
  } else {
    iprdamper = false;
  }
}

void md()
{

  // Solving equation of motion (a)
  integral_a[integral_type]();

  // Force and energy calculation
  potential(); if (mdmotion == 0) { return; }
  // External load
  loading();
  // Stretch
  ex=ex+dexdt*dt*1e12;
  ey=ey+deydt*dt*1e12;
  ez=ez+dezdt*dt*1e12;
  double strx=(1+ex)/(1+ex-dexdt*dt*1e12);
  double stry=(1+ey)/(1+ey-deydt*dt*1e12);
  double strz=(1+ez)/(1+ez-dezdt*dt*1e12);
  stretch(strx,stry,strz);

  if (repeat_lz) {
    if ((cell.hmat[2][2]*1e10 > repeat_lz_max)&&(dezdt>0)) { dezdt = -dezdt; }
    if ((cell.hmat[2][2]*1e10 < repeat_lz_min)&&(dezdt<0)) { dezdt = -dezdt; }
  }

  // Restrain cell change
  if (cell.pbcx == 0) {
    for (int i=0; i<3; i++) { cell.hmat[0][i] = cell.hmat_org[0][i]; } }
  if (cell.pbcy == 0) {
    for (int i=0; i<3; i++) { cell.hmat[1][i] = cell.hmat_org[1][i]; } }
  if (cell.pbcz == 0) {
    for (int i=0; i<3; i++) { cell.hmat[2][i] = cell.hmat_org[2][i]; } }

  // Solving equation of motion (b)
  integral_b[integral_type]();

  // Potential and kinetic energies
  atom.epotsum=0.0; atom.ekinsum=0.0;
  for (int i=1; i<=atom.natom; i++)
    {
      atom.epotsum=atom.epotsum+atom.epot[i];
      atom.ekinsum=atom.ekinsum+atom.wm[i]/2.0*(atom.vx[i]*atom.vx[i]+atom.vy[i]*atom.vy[i]+atom.vz[i]*atom.vz[i]);
    }

  // Velocity control
  if (ivscale == 1) { // Velocity scaling NVT
    vscale();
  } else if (irelax == 1) { 
    if (relax_algo == 0) {
      relax_gloc(); // GLOC relaxation
    } else {
      relax_fire();
    }
  }

  if (notrans) {
    double transx = 0.0; double transy = 0.0; double transz = 0.0;
    for (int i=1; i<=atom.natom; i++) {
      transx = transx + atom.vx[i]; transy = transy + atom.vy[i]; transz = transz + atom.vz[i];
    }
    transx = transx / atom.natom; transy = transy / atom.natom; transz = transz / atom.natom;
    for (int i=1; i<=atom.natom; i++) {
      atom.vx[i] = atom.vx[i] - transx; atom.vy[i] = atom.vy[i] - transy; atom.vz[i] = atom.vz[i] - transz;
    }
  }
  
  // For temperature monitoring
  if (atom.QC==0) {
    tempc = atom.Enkin()*2.0/3.0/(double)atom.natom/1.380662e-23;
  } else {
    tempc = atom.Enkin()*2.0/3.0/(double)atom.nrepatom/1.380662e-23;
  }

  // For output in "Status" area
  cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
  f_max=atom.Fmax()/ev*ang;
  dtm=dt*1.0e15;
  epotatom=atom.epotsum/ev/atom.natom;

  // For termination of MD
  bool recipe_prompt = false;
  if (itolfor) {
    if (f_max < tolfor) {
      mdmotion = 0; recipe_prompt = true; } }
  if (itolstep) {
    if (istep-istep0 >= tolstep) {
      mdmotion = 0; recipe_prompt = true;} }
  if ((irecipe)&&(recipe_prompt)) {
    recipe();
  }

  // Output config data
  /*
  if (confwrint > 0) {
    if (istep % confwrint == 0) {
      int num = istep / confwrint; 
      char filepath[80] = "CONFIG.SNAP"; char numc[10];
      if (num<10000) { snprintf(numc, sizeof(numc), "%04d", num);
      } else { snprintf(numc, sizeof(numc), "%d", num); }
      strcat(filepath,numc);
      writeconfig(filepath);
      if (autocap) { capture(); }
    } }
  */
}
  
