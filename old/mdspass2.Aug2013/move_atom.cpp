#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"
#include <GL/glui.h>

void potential();
void bookkeep();
void pot_initialize_all();
void set_atom_weight();
void set_atom_color();
int atom_number(char* at);
void deallocate_arrays();
void allocate_arrays();
extern GLuint objects;
extern GLfloat **color;
extern GLfloat yellow[];

void move_atom(int i, float frx, float fry, float frz)
{
  if ((i>0)&&(i<=atom.natom)) {
    atom.rx[i] = frx*ang;
    atom.ry[i] = fry*ang;
    atom.rz[i] = frz*ang;
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
  }
  pot_initialize_all();
  bookkeep();
  potential();
}

void move_atom(int i, const char *aasp, float frx, float fry, float frz)
{
  if ((i>0)&&(i<=atom.natom)) {
    atom.rx[i] = frx*ang;
    atom.ry[i] = fry*ang;
    atom.rz[i] = frz*ang;
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
    strcpy(atom.asp[i], aasp);
    atom.anum[i] = atom_number(atom.asp[i]);
  }
  set_atom_color();
  set_atom_weight();
  pot_initialize_all();
  bookkeep();
  potential();
}


