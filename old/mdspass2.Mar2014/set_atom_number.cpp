#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"
#include <GL/glui.h>

extern GLfloat **color, **color0;

int atom_number(char* at);

void set_atom_number()
{
  GLfloat yellow[] = {1.0, 1.0, 0.0, 1.0};
  GLfloat x[] = { 0.77, 1.0, 0.42 };
  GLfloat h[] = { 0.0, 0.95, 0.95 };
  GLfloat he[] = { 0.95, 0.95, 0.95 };
  GLfloat li[] = { 0.95, 0.95, 0.95 };
  GLfloat be[] = { 0.71, 0.71, 0.71 };
  GLfloat b[] = { 0.71, 0.71, 0.71 };
  GLfloat c[] = { 0.95, 0.95, 0.0 };
  GLfloat n[] = { 0.64, 0.80, 0.86 };
  GLfloat o[] = { 0.70, 0.0, 0.0 };
  GLfloat f[] = { 0.83, 0.95, 0.95 };
  GLfloat ne[] = { 0.95, 0.95, 0.95 };
  GLfloat na[] = { 0.0, 0.95, 0.95 };
  GLfloat mg[] = { 0.0, 0.95, 0.95 };
  GLfloat al[] = { 0.95, 0.0, 0.95 };
  GLfloat si[] = { 0.00, 0.95, 0.95 };
  GLfloat p[] = { 0.95, 0.95, 0.0 };
  GLfloat s[] = { 0.95, 0.95, 0.45 };
  GLfloat cl[] = { 0.71, 1.00, 0.00 };
  GLfloat ar[] = { 0.95, 0.95, 0.95 };
  GLfloat k[] = { 0.0, 0.95, 0.95 };
  GLfloat ca[] = { 0.0, 0.95, 0.95 };
  GLfloat sc[] = { 0.0, 0.95, 0.95 };
  GLfloat ti[] = { 0.60, 0.60, 0.60 };
  GLfloat v[] = { 0.64, 0.80, 0.86 };
  GLfloat cr[] = { 0.64, 0.80, 0.86 };
  GLfloat mn[] = { 0.64, 0.80, 0.86 };
  GLfloat fe[] = { 0.95, 0.00, 0.00 };
  GLfloat co[] = { 0.64, 0.80, 0.86 };
  GLfloat ni[] = { 0.64, 0.80, 0.86 };
  GLfloat cu[] = { 0.82, 0.45, 0.14 };
  GLfloat zn[] = { 0.95, 0.00, 0.95 };
  GLfloat ga[] = { 0.95, 0.00, 0.95 };
  GLfloat ge[] = { 0.95, 0.00, 0.95 };
  GLfloat as[] = { 0.95, 0.95, 0.00 };
  GLfloat se[] = { 0.95, 0.95, 0.00 };
  GLfloat br[] = { 0.95, 0.00, 0.00 };
  GLfloat kr[] = { 0.75, 0.75, 0.75 };
  GLfloat rb[] = { 0.00, 0.95, 0.95 };
  GLfloat sr[] = { 0.00, 0.95, 0.95 };
  GLfloat y[] = { 0.25, 0.95, 0.25 };
  //GLfloat y[]  = { 0.75, 0.75, 0.75 };
  GLfloat zr[] = { 0.75, 0.75, 0.75 };
  GLfloat nb[] = { 0.75, 0.75, 0.75 };
  GLfloat mo[] = { 0.75, 0.75, 0.75 };
  GLfloat tc[] = { 0.75, 0.75, 0.75 };
  GLfloat ru[] = { 0.75, 0.75, 0.75 };
  GLfloat rh[] = { 0.75, 0.75, 0.75 };
  GLfloat pd[] = { 1.00, 1.00, 1.00 };
  GLfloat ag[] = { 1.00, 1.00, 1.00 };
  GLfloat cd[] = { 0.75, 0.75, 0.75 };
  GLfloat in[] = { 0.75, 0.75, 0.75 };
  GLfloat sn[] = { 0.75, 0.75, 0.75 };
  GLfloat sb[] = { 0.75, 0.75, 0.75 };
  GLfloat te[] = { 0.75, 0.75, 0.75 };
  GLfloat I[] = { 0.75, 0.75, 0.75 };
  GLfloat xe[] = { 0.75, 0.75, 0.75 };
  GLfloat cs[] = { 0.00, 0.95, 0.95 };
  GLfloat ba[] = { 0.00, 0.95, 0.95 };
  GLfloat la[] = { 0.75, 0.75, 0.75 }; //57
  //ce,pr,nd,pm,sm,eu,gd,tb,dy,ho,er,tm,yb,lu,hf,ta,w,re,os,ir
  GLfloat pt[] = { 0.75, 0.75, 0.75 }; //78
  GLfloat au[] = { 1.00, 0.85, 0.00 };
  GLfloat hg[] = { 0.75, 0.75, 0.75 }; //80
  //tl,pb,bi,po,at,rn,fr,ra,ac,th,pa,u,np,pu,am,cm,bk
  GLfloat cf[] = { 0.75, 0.75, 0.75 }; //98
  GLfloat es[] = { 0.93, 0.88, 0.80 };
  GLfloat fm[] = { 0.93, 0.88, 0.80 };
  for (int i=1; i<=atom.natom*3; i++) { memcpy(color[i],x,sizeof(GLfloat)*4); }
  for (int i=1; i<=atom.natom; i++) {
    if (atom_number(atom.asp[i])==1) memcpy(color[i],h,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==2) memcpy(color[i],he,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==3) memcpy(color[i],li,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==4) memcpy(color[i],be,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==5) memcpy(color[i],b,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==6) memcpy(color[i],c,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==7) memcpy(color[i],n,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==8) memcpy(color[i],o,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==9) memcpy(color[i],f,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==10) memcpy(color[i],ne,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==11) memcpy(color[i],na,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==12) memcpy(color[i],mg,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==13) memcpy(color[i],al,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==14) memcpy(color[i],si,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==15) memcpy(color[i],p,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==16) memcpy(color[i],s,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==17) memcpy(color[i],cl,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==18) memcpy(color[i],ar,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==19) memcpy(color[i],k,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==20) memcpy(color[i],ca,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==21) memcpy(color[i],sc,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==22) memcpy(color[i],ti,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==23) memcpy(color[i],v,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==24) memcpy(color[i],cr,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==25) memcpy(color[i],mn,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==26) memcpy(color[i],fe,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==27) memcpy(color[i],co,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==28) memcpy(color[i],ni,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==29) memcpy(color[i],cu,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==30) memcpy(color[i],zn,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==31) memcpy(color[i],ga,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==32) memcpy(color[i],ge,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==33) memcpy(color[i],as,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==34) memcpy(color[i],se,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==35) memcpy(color[i],br,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==36) memcpy(color[i],kr,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==37) memcpy(color[i],rb,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==38) memcpy(color[i],sr,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==39) memcpy(color[i],y,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==40) memcpy(color[i],zr,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==41) memcpy(color[i],nb,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==42) memcpy(color[i],mo,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==43) memcpy(color[i],tc,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==44) memcpy(color[i],ru,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==45) memcpy(color[i],rh,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==46) memcpy(color[i],pd,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==47) memcpy(color[i],ag,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==48) memcpy(color[i],cd,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==49) memcpy(color[i],in,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==50) memcpy(color[i],sn,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==51) memcpy(color[i],sb,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==52) memcpy(color[i],te,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==53) memcpy(color[i],I,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==54) memcpy(color[i],xe,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==55) memcpy(color[i],cs,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==56) memcpy(color[i],ba,sizeof(GLfloat)*4);
    if ((atom_number(atom.asp[i])>=57)&&(atom_number(atom.asp[i])<=78))
      memcpy(color[i],la,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==79) memcpy(color[i],au,sizeof(GLfloat)*4);
    if ((atom_number(atom.asp[i])>=80)&&(atom_number(atom.asp[i])<=98))
      memcpy(color[i],hg,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==99) memcpy(color[i],es,sizeof(GLfloat)*4);
    if (atom_number(atom.asp[i])==100) memcpy(color[i],fm,sizeof(GLfloat)*4);
  }
  for (int i=1; i<=atom.natom*3; i++) { memcpy(color0[i],color[i],sizeof(GLfloat)*4); }
}
