#include<math.h>
#include<string.h>
#define NMAX 2500
#define NBOOK 1000
#define ang 1.0e-10
#define ev  1.6021892e-19
#define GPa 1.0e9
#define MPa 1.0e6
#define SELECTIONS 10
#define MAXPOTTYPE 8
#define MAXENSTYPE 9
#define MAXNEIGHBOR 100
#define NINTEGRAL 10
#define MAXMODE 30
#define MAXKP 5

#define FLOAT double
#define FLOAT2 double

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include "pair_function.h" // Kubo 20140224

class Atom{
 public:
  Atom();
  int natom, nrepatom;
  char **asp; //Atom species
  int *anum; //Atomic number
  char potential_func[30];
  char potential_arg[30];
  double *wm;
  double *rx, *ry, *rz, *fx, *fy, *fz;
  double *rx_org, *ry_org, *rz_org, *rx_p, *ry_p, *rz_p;
  double *vx, *vy, *vz;
  double *ax, *ay, *az, *bx, *by, *bz, *cx, *cy, *cz; // for p-c integral
  double *fx_l, *fy_l, *fz_l; //loading
  double epotsum, ekinsum;
  double *epot, *epot_p;
  FLOAT *rx_float, *ry_float, *rz_float, *fx_float, *fy_float, *fz_float;
  FLOAT2 *epot_float;
  bool *mfx, *mfy, *mfz; // To fix atom motion
  bool *lock; // To lock atom (for gloc-accel)
  int nelem; // number of elements
  int *repatom; // whether an atom is repatom or not (1=yes, 0=no)
  int *elem_id; // element number to which an atom belongs
  int **elem_v, **elem_v_rep;
  int nwall; // number of wall-triangle (for CNT)
  int **wall_v; // vertices of wall (for CNT)  ==> wall_v[nwall][3]
  double *wall_area; // area of wall (for CNT) ==> wall_area[nwall]
  int *wall_vrev; // reverse normal vector or not
  double **wall_nvec; // normal vector of wall (for CNT) ==> wall_nvec[nwall][3]
  int **neighbor; // list of neiboring atoms
  int *nneighbor; // number of neiboring atoms
  double Enkin();
  double Dist2(int i, int j);
  double Dist2(int i, int j, int ix, int iy, int iz);
  double Dist2Closest(int i, int j, int &ix, int &iy, int &iz);
  double Dx(int i, int j, int ix, int iy, int iz);
  double Dy(int i, int j, int ix, int iy, int iz);
  double Dz(int i, int j, int ix, int iy, int iz);
  double Angle(int i, int j, int k, int ix0, int iy0, int iz0, int ix1, int iy1, int iz1);
  double Fmax();
  int QC;
  int instcenter;
  double **evecx, **evecy, **evecz;
  double eigval[MAXMODE];
  double ***satom;
  int getAtomNumber(char* atom){
    if ((strcmp(atom,"C")==0)||(strcmp(atom,"c")==0)) { return 6;
    } else if ((strcmp(atom,"Cu")==0)||(strcmp(atom,"cu")==0)) { return 29;
    } else if ((strcmp(atom,"Al")==0)||(strcmp(atom,"al")==0)) { return 13;
    } else if ((strcmp(atom,"Si")==0)||(strcmp(atom,"si")==0)) { return 14;
    } else if ((strcmp(atom,"Ni")==0)||(strcmp(atom,"ni")==0)) { return 22;
    } else if ((strcmp(atom,"Fe")==0)||(strcmp(atom,"fe")==0)) { return 26;
    } else if ((strcmp(atom,"Ti")==0)||(strcmp(atom,"ti")==0)) { return 28;
    } else if ((strcmp(atom,"Ge")==0)||(strcmp(atom,"ge")==0)) { return 32;
    } else if ((strcmp(atom,"Ag")==0)||(strcmp(atom,"ag")==0)) { return 47;
    } else if ((strcmp(atom,"Sn")==0)||(strcmp(atom,"sn")==0)) { return 50;
    } else if ((strcmp(atom,"Au")==0)||(strcmp(atom,"au")==0)) { return 79;
    } else if ((strcmp(atom,"Pd")==0)||(strcmp(atom,"pd")==0)) { return 46;
    } else if ((strcmp(atom,"Pt")==0)||(strcmp(atom,"pt")==0)) { return 78;
    } else { return 1; }
  }
 private:
  double dist2;
};

class Book{
 public:
  Book();
  int *alistnum;
  int ***alist;
  double frc, frc2;
  int nbk;
  int algo;
  int natom, nbook, nbook_new;
  bool alloc;
};

class Cell{
 public:
  Cell();
  double hmat[3][3], hinmat[3][3], hvmat[3][3];
  double hamat[3][3], hbmat[3][3], hcmat[3][3];
  double hmat_org[3][3];
  double len[3];
  int pbcx, pbcy, pbcz;
  double alat;
  double volume;
  double virx, viry, virz;
  double sgmmat[3][3], dmat[3][3], sgmmat_set[3][3];
  double sgmmat_p[3][3], hmat_p[3][3];
  bool relax_static_initial;
  double ecmat[6][6];
  double Getvolume();
  void Reset();
  void Setlen();
  double Dstress();
  double ww, prmass, damper_param;
  bool fix[3][3];
};

class Tersoff{
 public:
  Tersoff() {
    initialize = true;
  };
  double **zmat, **b;
  double terrr, terss, teraa, terbb, terlambda, termu, terbeta;
  double tern, terc, terd, terh, terw;
  int term;
  double terc2, terd2, termum;
  bool initialize;
};

class GEAM{
 public:
  GEAM();
  double gre[16], gfe[16], grhoe[16], grhos[16], galp[16], gbet[16], gaa[16], gbb[16], gkai[16],
    glam[16], gffn0[16], gffn1[16], gffn2[16], gffn3[16], gff0[16], gff1[16], gff2[16], gff3[16],
    geta[16], gffe[16];
  int *sp;
  double *rhob;
  bool initialize;
};

class Eammis{
 public:
  Eammis() {
    initialize = true;
    mesh = 10000;
    //mesh = 1000000;
  };
  double *pair, *paird, *den, *dend, *embed, *embedd;
  double pairmin, pairmax, denmin, denmax, embedmin, embedmax;
  double *pairx, *pairy, **pairc, *pairxn;
  double *denx, *deny, **denc, *denxn;
  double *embedx, *embedy, **embedc, *embedxn;
  double *rhob;
  int mesh;
  bool initialize;
};

class Adp{
 public:
  Adp() { // Constructor
    initialize = true; cut = 7.0;
#if defined __linux__ || defined __APPLE__
    strcpy(fname,"pot/ADP_NdFeB.pot");
#else
    strcpy(fname,"pot\\ADP_NdFeB.pot");
#endif
  };
  // Potential parameters should come here
  double cut;
  bool initialize;
  double *rho, *mu_x, *mu_y, *mu_z, *nu;
  double *lambda_xx, *lambda_yy, *lambda_zz;
  double *lambda_xy, *lambda_yz, *lambda_zx;
  //double *param_pair, *param_dp, *param_qp;
  //double *param_den, *param_eb, *gradF;
  double **param_pair, **param_dp, **param_qp;
  double **param_den, **param_eb, *gradF;
  int ntype, nptype;
  int *type, **ptype, *typen;
  char fname[60];
// Kubo 20140224 -------------//
  PairFunction **Phi;
  PairFunction **Rho;
  PairFunction **F;
  PairFunction **U;
  PairFunction **W;
// End -----------------------//
};

class Dipole{
 public:
  Dipole() {
    dp_eps = 14.40; dp_cut = 8.0; dp_tol = 1.e-7; dp_mix = 0.2;
    // ntype = 3; nptype = 3;
    initialize = true;
    shift = true;
//    shift = false; // Kubo
#if defined __linux__ || defined __APPLE__
    strcpy(fname,"pot/Dipole_YSZ.pot");
#else
    strcpy(fname,"pot\\Dipole_YSZ.pot");
#endif
  };
  double *dp_alpha, *dp_b, *dp_c, *r_cut;
  double *E_statx, *E_staty, *E_statz;
  double *E_indx, *E_indy, *E_indz;
  double *E_oldx, *E_oldy, *E_oldz;
  double *E_totx, *E_toty, *E_totz;
  double *p_srx, *p_sry, *p_srz;
  double *p_indx, *p_indy, *p_indz;
  double dp_eps, dp_cut, dp_tol, dp_mix, ew_rcut;
  double *ratio, *charge, last_charge, dp_kappa;
  double *buck_a, *buck_s, *buck_c;
  double *buck_vrc, *buck_vprc;
  int sw_kappa;
  bool initialize, shift;
  int ntype, nptype;
  int *type, **ptype, *typen;
  char fname[60];
  double val1,val2,val3;
};

class Bre{
 public:
  Bre() {
    //    std::cout<< "Bre constructor was called" << std::endl;
    initialize = true;
    npmax = 7000; nlmax = 1000000, ntab = 10000;
    npm1 = npmax + 1; ntypes = 8; nnma = 5000*40;
    ikl = 2; nmabig = nnma*5; ver = 3.00;
    xqm = 3.70; att = 3.20;
    ljmin = 2.0; ljmax = 10.0; ljmesh = 10000;
    //    std::cout<<nmabig<<std::endl;
    };
  int npmax, nlmax, ntab, npm1, ntypes, nnma, ikl,
    nmabig;
  double ver;
  int igh[26], in3[65][4], ndihed;
  double *bww, *dww, *rcor, *cor1, *cor2, *cor3;
  int *rep1, *rep2, *rep3;
  double xh[3][11][11], xh1[3][11][11], xh2[3][11][11];
  int *list, *nabors, *lcheck;
  double ad[5][5], axl[5][5], bd[5][5], bxl[5][5],
    cd[5][5], cxl[5][5], dd[5][5], dxl[5][5], ed[5][5],
    xn1[3], rb1[5][5], rb2[5][5], pid[5][5], rmax[5][5],
    rlist[5][5],spgc[7][6],spgh[7][4], xq, att, xqm, pq,
    xtn2[5], xtn1[5], adb[5], cdb[5], cd2[5], ddb[5],
    ddb2[5], hdb[5], chi[5][5], xtm[5], xtl[5];
  double epslj[5][5], siglj[5][5], epsts[5][5], bmin[5][5], bmax[5][5]; //AIREBO
  int igc[26];
  double xxdb, xdb[5][5][5], reg[5][5][5];
  double clmn[4][11][11][11][65], clm[3][11][11][17],
    tlmn[11][11][11][65], pidt;
  int in2[17][3];
  double *exx1, *dexx1, *exx2;
  int *ivct2b, *jvct2b;
  double *xhc1, *xhc2;
  double ***atable, ***datable, ***rtable, ***drtable,
    ddtab[5][5], ***tabfc, ***tabdfc;
  bool initialize;

  double **eps, **sig, *xmass;
  int *noa;
  double rhh, rch, sigma, epsi, rll;
  double rc, rlis;
  int *ktype, kt[101], kt2[101], np;
  double *r01, *r02, *r03, *rnp1, *rnp2, *rnp3;
  double *rpp1, *rpp2, *rpp3;
  double tote;
  int kend;
  //AIREBO
  double *lj, *ljd, ljmin, ljmax;
  int ljmesh;
};

class Cond{
 public:
  Cond();
};

class Fire{
 public:
  Fire();
  double sp,fnorm,vnorm,ffinc,ffdec,ffalph,fire_alph,fire_alph_ini,ffdtmax;
  int ifire,nfmin;
};

class Cg{
 public:
  Cg();
  double lamtr, lam_init;
  double **sgx, **sgy, **sgz, **shx, **shy, **shz;
  double *rstx, *rsty, *rstz;
  double pot0, pot1, grad0, grad1, normsh;
  int step;
};

extern int istep;
extern int mdmotion;
extern double rcut, rcut2;
extern float rcut_f, frcmar, frc_f;
extern float temp_set;
extern float tempc,cellx,celly,cellz,f_max,epotatom;
extern float cell1x,cell1y,cell1z,cell2x,cell2y,cell2z,cell3x,cell3y,cell3z;
extern float slice_1min,slice_1max,slice_2min,slice_2max,slice_3min,slice_3max;
extern int ifix_atoms;
extern float strs_xx,strs_yy,strs_zz,strs_xy,strs_yz,strs_zx;
extern int ensemble;
extern int config_type;
extern int irepx, irepy, irepz, icntm, icntn;
extern float alat, cscnt, rotz, shiftz;
extern int imerge;
extern char current_config_name[100], current_qcelement_name[100];
extern const char  *potstring_list[];
extern int ipottype;
extern int notrans, incell;
extern FILE *energyfile, *stressfile, *cellfile, *ssnorfile, *ecnorfile, *ecallfile;
extern double dt;
extern int repeat_lz;
extern float dexdt, deydt, dezdt, repeat_lz_min, repeat_lz_max;
extern int yzplane_punch;
extern float yzplane_punch_d, yzplane_punch_dd, yzplane_punch_ftot;
//extern float ex;
//extern float dexdt;
extern void (*integral_a[NINTEGRAL])(), (*integral_b[NINTEGRAL])();
extern int integral_type;
extern float prmass_scale;
extern int initial_velocity;

extern Atom atom;
extern Cell cell;
extern Book book;
extern Tersoff tersoff;
extern GEAM geam;
extern Eammis eammis;
extern Dipole dipole;
extern Cond cond;
extern Bre bre;
extern Fire fire;
extern Adp adp;
extern Cg cg;

#ifdef CUDA
extern  int *arrayDint;
extern  int *arrayDrepatom;
extern  double *arrayDrx, *arrayDry, *arrayDrz, *arrayDfx, *arrayDfy, *arrayDfz;
extern  double *arrayDepot, *arrayDhmat, arrayHhmat[9];
extern  int *arrayDalistnum, *arrayDalist;
extern  int iarray[(NMAX+1)*(NBOOK+1)*4];
#endif
