#include <fstream>
#include <iostream>
#include "myheader.h"
#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

void capture();
void writeconfig(const char* fname);
extern int confwrint, autocap;

/*
void writedata()
{
  std::ofstream foutene( "energy.d", std::ios::out | std::ios::app );
  foutene << istep << " " << atom.epotsum << " " << atom.ekinsum << std::endl;
  
}
*/
//void writedata( FILE *fp )
void writedata()
{
  fprintf(energyfile, "%d %22.17e %22.17e %22.17e\n",
	  istep, atom.epotsum/ev/atom.natom, atom.ekinsum/ev/atom.natom,
	  (atom.epotsum+atom.ekinsum)/ev/atom.natom);
  fprintf(stressfile, "%d %10.5e %10.5e %10.5e  %10.5e %10.5e %10.5e\n",
	  istep, cell.sgmmat[0][0]*1e-6, cell.sgmmat[1][1]*1e-6, cell.sgmmat[2][2]*1e-6,
	  cell.sgmmat[0][1]*1e-6, cell.sgmmat[1][2]*1e-6, cell.sgmmat[2][0]*1e-6);
  fprintf(ssnorfile, "%10.5e %10.5e %10.5e  %10.5e %10.5e %10.5e\n",
	  cell.hmat[0][0]*1e10, cell.hmat[1][1]*1e10, cell.hmat[2][2]*1e10,
	  cell.sgmmat[0][0]*1e-6, cell.sgmmat[1][1]*1e-6, cell.sgmmat[2][2]*1e-6);
  fprintf(cellfile, "%d %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f\n",
	  istep,
	  cell.hmat[0][0]*1e10, cell.hmat[1][0]*1e10, cell.hmat[2][0]*1e10,
	  cell.hmat[0][1]*1e10, cell.hmat[1][1]*1e10, cell.hmat[2][1]*1e10,
	  cell.hmat[0][2]*1e10, cell.hmat[1][2]*1e10, cell.hmat[2][2]*1e10);
  fflush(energyfile); fflush(stressfile); fflush(cellfile); fflush(ssnorfile);

  // Output config data
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

  /*
  for (int i=1;i<=atom.natom;i++) {
    printf("%d %22.17e\n",i,atom.epot[i]);
  }
  */
}
/* void writedata_initialize()
{
  std::ofstream foutene( "energy.d", std::ios::trunc );
} */
 //int writedata_initialize(FILE *fp)
void writedata_initialize()
{
  if (energyfile != NULL) { fclose(energyfile); }
  energyfile = fopen("energy.d","w");
  fprintf(energyfile, "# step, Epot (eV, per atom), Ekin, Epot+Ekin\n");
  if (stressfile != NULL) { fclose(stressfile); }
  stressfile = fopen("stress.d","w");
  fprintf(stressfile, "# step, Stress xx (MPa), yy, zz, xy, yz, zx\n");
  if (ssnorfile != NULL) { fclose(ssnorfile); }
  ssnorfile = fopen("ssnormal.d","w");
  fprintf(ssnorfile, "# cell 11 22 33 Stress xx (MPa), yy, zz\n");
  if (cellfile != NULL) { fclose(cellfile); }
  cellfile = fopen("cell.d","w");
  fprintf(cellfile, "# step, Cell matrix 11 (A), 21, 31,  12, 22, 32,  13, 23, 33\n");
}
