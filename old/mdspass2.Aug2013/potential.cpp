#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

void e_force_morse();
void e_force_geam_am();
void e_force_brenner();
void e_force_brenner(int mode);
void e_force_tersoff();
void e_force_mishin();
void e_force_dipole();
void e_force_adp();
void calcstresstensor();
void bookkeep();
void get_first_arg(std::string &line, std::string &arg1);
int count_arg_number(std::string line);

extern std::string potfile;

void potential_set (int control)
{
  printf("Potential: No. %d, %s\n",ipottype,potstring_list[ipottype]);
  //strcpy(atom.potential_func, potstring_list[ipottype]);
  std::string line = std::string(potstring_list[ipottype]);
  std::string arg1, arg2; int narg = count_arg_number(line);
  if (narg==1) { get_first_arg(line,arg1);
    strcpy(atom.potential_func, arg1.c_str());
    strcpy(atom.potential_arg, "");
  } else { get_first_arg(line,arg1); get_first_arg(line,arg2);
    strcpy(atom.potential_func, arg1.c_str());
    strcpy(atom.potential_arg,  arg2.c_str()); }
  book.algo = 1;
  if (strcmp(atom.potential_func,"Tersoff")==0) { book.algo = 2; }
  if (strcmp(atom.potential_func,"Brenner")==0) { book.algo = 3; }
  if (strcmp(atom.potential_func,"AIREBO" )==0) { book.algo = 3; }
}

void potential()
{
 CALC:
  if (strcmp(atom.potential_func, "Morse") == 0) {
    e_force_morse();
  } else if (strcmp(atom.potential_func, "GEAM") == 0) {
    e_force_geam_am();
  } else if (strcmp(atom.potential_func, "Brenner") == 0) {
    e_force_brenner();
  } else if (strcmp(atom.potential_func, "AIREBO") == 0) {
    e_force_brenner(1);
  } else if (strcmp(atom.potential_func, "Tersoff") == 0) {
    e_force_tersoff();
  } else if ((strcmp(atom.potential_func, "EAM") == 0)&&
	     (strcmp(atom.potential_arg, "Mishin") == 0)) {
    e_force_mishin();
  } else if (strcmp(atom.potential_func, "Dipole") == 0) {
    e_force_dipole();
  } else if (strcmp(atom.potential_func, "ADP") == 0) {
    e_force_adp();
  } else {
    printf("##### Invalid potential for this system #####\n");
    mdmotion = 0;
    return;
  }

  // Set rcut_f
  if (log10(rcut)<-8) { rcut_f = rcut*1e10; } // Some routines use rcut in [m]
  else { rcut_f = rcut; }

  // Is frc enough?
  //printf("%f %f\n",rcut_f, book.frc/ang);
  if (rcut_f > book.frc / ang) {
    book.frc = ((double)rcut_f + frcmar) * ang;
    frc_f = book.frc / ang;
    bookkeep();
    goto CALC;
  }
  // For output in "Stresss" area
  calcstresstensor();
  strs_xx = cell.sgmmat[0][0]*1e-6;
  strs_yy = cell.sgmmat[1][1]*1e-6;
  strs_zz = cell.sgmmat[2][2]*1e-6;
  strs_xy = cell.sgmmat[0][1]*1e-6;
  strs_yz = cell.sgmmat[1][2]*1e-6;
  strs_zx = cell.sgmmat[2][0]*1e-6;
}

void set_potfile()
{
  if (strcmp(atom.potential_func, "Dipole") == 0) {
    printf("Dipole parameter file is set to %s\n",potfile.c_str());
    strcpy(dipole.fname,potfile.c_str());
    dipole.initialize = true;
  }
}

void pot_initialize_all()
{
  bre.initialize = true;
  tersoff.initialize = true;
  eammis.initialize = true;
  geam.initialize = true;
  adp.initialize = true;
  dipole.initialize = true;
}
