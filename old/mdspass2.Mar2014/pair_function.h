#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

class PairFunction{ //========================================//


 protected: //==================================================

//////// Member Variables //////////////////////////////////////
  int     NumOfParam;   // Number of Parameters in Function.
  double  CutoffRadius; // Cutoff Radius for Interaction. Nonpositive Value kills Cutoff.
  double *ParamTable;   // Potential Paramater Table
  void  (*Func)(double, PairFunction*, double&, double&);
  string  ID;           // Identifier for Error Message (Unimplemented...).

  bool FlgWarningForNumericalDerivative;

 private: //====================================================

  // Initializer ---------------
  // Initialize Members by naught.
  void Initialize(void){
    NumOfParam = 0;
    CutoffRadius = 0.0;
    ParamTable = NULL;
    Func = NULL;
    ID ="";

    FlgWarningForNumericalDerivative = true;
  }

  // Set NumOfParam, Func, ParamTable.
  void Initialize(string funcform, double *param){
    Initialize();

    Func = GetFunc(funcform);
    NumOfParam = GetNumOfParam(funcform);
    if(NumOfParam>0){
      ParamTable = new double[NumOfParam];
      for(int i=0;i<NumOfParam;i++){ParamTable[i] = param[i];}
    }
  }

  // Set Above + CutoffRadius.
  void Initialize(string funcform, double *param, double rc){
    Initialize(funcform,param);
    CutoffRadius = rc;
  }

  // Set Above + ID.
  void Initialize(string funcform, double *param,double rc,string tag){
    Initialize(funcform,param,rc);
    ID = tag;
  }

  // Numerical Derivative
  double NumericalDerivative(double rr,double dr){
    if(FlgWarningForNumericalDerivative){
      cout<<"Numerical Derivative Applied";
      if(ID==""){cout<<"."<<endl;}
      else{cout<<" for "<< ID << "."<<endl;}
      FlgWarningForNumericalDerivative = false;
    }

    double fx0,fx1,fx2;
    double dfdx0,dfdx1,dfdx2;

    (*Func)(rr,   this,fx0,dfdx0);
    (*Func)(rr+dr,this,fx1,dfdx1);
    (*Func)(rr-dr,this,fx2,dfdx2);

    dfdx0 = (fx1-fx2)/(2.0*dr);
    return(dfdx0);
  }

 public: //=====================================================
  // Constructor ---------------
  PairFunction(void){
    Initialize();
  }

  PairFunction(string funcform, double *param){
    Initialize(funcform,param);
  }

  PairFunction(string funcform, double *param, double rc){
    Initialize(funcform,param,rc);
  }

  PairFunction(string funcform, double *param, double rc,string s){
    Initialize(funcform,param,rc,s);
  }


  // Destructor ----------------
 ~PairFunction(){
    if(ParamTable!=NULL){
      delete [] ParamTable;
      ParamTable = NULL;
    }
  }

//////// Common Functions //////////////////////////////////////
  // Calculate Function
  void Calc(double rr, double &val, double &grad){
    if((CutoffRadius>0.0)   // Cutoff Working?
     &&(rr>CutoffRadius)){  // Out of Cutoff?
      val  = 0.0;
      grad = 0.0;
    }else{
      (*Func)(rr,this,val,grad);
//      grad = NumericalDerivative(rr,0.001); // For Debug;
    }
  }

  // Plotter for Debugging
  void Plot(
    double min,  // Minmum  Value to Plot
    double max,  // Maximum Value to Plot
    int nmesh,   // Number of Mesh Data
    string fname // File Name to Output
  ){
    double x,fx,dfdx,dfdx_nmr;
    double dx = fabs(max-min)/nmesh;
    ofstream ofs(fname.c_str());

    if(min>max){min = max;}

    ofs << "# r, f(r), df/dr(r)_analytic, df/dr(r)_numerical" << endl;

    bool flgtmp = FlgWarningForNumericalDerivative;
    FlgWarningForNumericalDerivative = false;

    for(int i=0;i<=nmesh;i++){
      x = min + dx*i;

      Calc(x,fx,dfdx);
      dfdx_nmr = NumericalDerivative(x,0.001);

      ofs << setiosflags(ios::scientific)
          << setw(16) << setprecision(6)  << x;
      ofs << setiosflags(ios::scientific)
          << setw(22) << setprecision(10) << fx;
      ofs << setiosflags(ios::scientific)
          << setw(22) << setprecision(10) << dfdx;
      ofs << setiosflags(ios::scientific)
          << setw(22) << setprecision(10) << dfdx_nmr;
      ofs << endl;
    }

    FlgWarningForNumericalDerivative = flgtmp;
  }

//////// Reference to Members //////////////////////////////////
  int GetNumOfParam(void){
    return(NumOfParam);
  }

  double GetParam(int id){
    if(id<NumOfParam){
      return(ParamTable[id]);
    }else{
      cout << "WARNING: Mismatch in Number of Parameters." << endl;
      return(0.0);
    }
  }





////////////////////////////////////////////////////////////////X
////                                                        ////X
////          EDIT BELOW TO ADD NEW FUNCTION FORMS          ////X
////                                                        ////X
////////////////////////////////////////////////////////////////X

  int GetNumOfParam(string funcform){
    if(funcform=="const")      {return(1);}
    if(funcform=="csw2")       {return(4);}
    if(funcform=="csw2_sc")    {return(5);}
    if(funcform=="eopp_exp")   {return(6);}
    if(funcform=="eopp_exp_sc"){return(7);}
    if(funcform=="exp_plus")   {return(3);}
    if(funcform=="exp_plus_sc"){return(4);}
    if(funcform=="meopp")      {return(7);}
    if(funcform=="meopp_sc")   {return(8);}
    if(funcform=="mishin")     {return(6);}
    if(funcform=="mishin_sc")  {return(7);}
    if(funcform=="ms")         {return(3);}
    if(funcform=="poly_5")     {return(5);}
    return(0);
  }

 protected: //==================================================
  void (*GetFunc(string funcform))(double, PairFunction*, double&, double&){
    if(funcform=="const")      {return(Const);}
    if(funcform=="csw2")       {return(Csw2);}
    if(funcform=="csw2_sc")    {return(Csw2_sc);}
    if(funcform=="eopp_exp")   {return(Eopp_Exp);}
    if(funcform=="eopp_exp_sc"){return(Eopp_Exp_sc);}
    if(funcform=="exp_plus")   {return(Exp_Plus);}
    if(funcform=="exp_plus_sc"){return(Exp_Plus_sc);}
    if(funcform=="meopp")      {return(Meopp);}
    if(funcform=="meopp_sc")   {return(Meopp_sc);}
    if(funcform=="mishin")     {return(Mishin);}
    if(funcform=="mishin_sc")  {return(Mishin_sc);}
    if(funcform=="ms")         {return(MS);}
    if(funcform=="poly_5")     {return(Poly_5);}
    cout << "Unknown Function ("<< funcform <<") Type Found." << endl;
    cout << "Zero Function, f(x)=0, Applied for First Aid." << endl;
    return(Zero);
  }

// PSI Function for Cutoff -------------------------------------
  void Psi(
        double rr,      //  IN: Distance (Angst.)
        double h,       //  IN: Smoothing Parameter (Angst.)
        double &val,    // OUT: Value (-)
        double &grad    // OUT: Derivative (1/Angst.)
  ){
    double x  = (rr-CutoffRadius)/h;
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x2*x2;
    if(rr<CutoffRadius){
      val  = x4/(1.0+x4);
      grad = 4.0*x3/(1.0+x4)/(1.0+x4)/h;
    }else{
      val  = 0.0;
      grad = 0.0;
    }
  }


// Template ----------------------------------------------------
  static void Func_Template(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    val  = pf->GetParam(0)*0.0;
    grad = pf->GetParam(0)*0.0;
  }
  static void Func_Template_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(0);
    Func_Template(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// Zero --------------------------------------------------------
  static void Zero(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    val  = 0.0;
    grad = 0.0;
  }


// Const -------------------------------------------------------
  static void Const(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    val  = pf->GetParam(0);
    grad = 0.0;
  }


// CSW2 --------------------------------------------------------
  static void Csw2(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double p = pow(rr, pf->GetParam(3));
    val  = (1.0+pf->GetParam(0)*cos(pf->GetParam(1)*rr+pf->GetParam(2)))/p;
    grad = (-pf->GetParam(0)*pf->GetParam(1)*sin(pf->GetParam(1)*rr+pf->GetParam(2)))/p
         - pf->GetParam(3)*val/rr;
  }
  static void Csw2_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(pf->GetNumOfParam()-1);
    Csw2(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// Eopp_Exp ----------------------------------------------------
  static void Eopp_Exp(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double tmp[5];
    tmp[0] = pf->GetParam(0)*exp(-pf->GetParam(1)*rr);
    tmp[1] = pf->GetParam(2)/pow(rr,pf->GetParam(3));
    tmp[2] = pf->GetParam(4)*rr + pf->GetParam(5);
    tmp[3] = cos(tmp[2]);
    tmp[4] = sin(tmp[2]);

    val  = tmp[0] + tmp[1]*tmp[3];
    grad = - pf->GetParam(1)*tmp[0]
           - pf->GetParam(3)*tmp[1]*tmp[3]/rr
           - pf->GetParam(4)*tmp[1]*tmp[4];
  }
  static void Eopp_Exp_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(pf->GetNumOfParam()-1);
    Eopp_Exp(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// Exp_Plus ----------------------------------------------------
  static void Exp_Plus(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double tmp = pf->GetParam(0)*exp(-pf->GetParam(1)*rr);
    val  = tmp + pf->GetParam(2);
    grad = -pf->GetParam(1) * tmp;
  }
  static void Exp_Plus_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(pf->GetNumOfParam()-1);
    Exp_Plus(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// Meopp -------------------------------------------------------
  static void Meopp(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double tmp[6];

    tmp[0] = rr - pf->GetParam(6);
    tmp[1] = pf->GetParam(0)/pow(tmp[0],pf->GetParam(1));
    tmp[2] = pf->GetParam(2)/pow(rr,pf->GetParam(3));
    tmp[3] = pf->GetParam(4)*rr + pf->GetParam(5);
    tmp[4] = cos(tmp[3]);
    tmp[5] = sin(tmp[3]);

    val  = tmp[1] + tmp[2]*tmp[4];
//    grad = - pf->GetParam(2)*tmp[1]/tmp[0]
    grad = - pf->GetParam(1)*tmp[1]/tmp[0] // Kubo 20140318
           - pf->GetParam(3)*tmp[2]*tmp[4]/rr
           - pf->GetParam(4)*tmp[2]*tmp[5];

  }
  static void Meopp_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(pf->GetNumOfParam()-1);
    Meopp(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// Mishin ------------------------------------------------------
  static void Mishin(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double z  = rr - pf->GetParam(3);
//    double e  = exp(-pf->GetParam(5)*z);
    double e  = exp(-pf->GetParam(5)*rr); // Defition in POTFIT
    double p  = pow(z,pf->GetParam(4));
    double pp = pf->GetParam(4)*p/z;
    double ep = -pf->GetParam(5)*e;
    val  = pf->GetParam(0)*  p*e*(1.0+pf->GetParam(1)*e) + pf->GetParam(2);
    grad = pf->GetParam(0)*(pp*e*(1.0+pf->GetParam(1)*e)
  			     + p*ep*(1.0+pf->GetParam(1)* e)
  			     + p* e*(    pf->GetParam(1)*ep));
  }
  static void Mishin_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(pf->GetNumOfParam()-1);
    Mishin(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// MS (Morse Stretch) ------------------------------------------
  static void MS(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){

    double de,a,r0;
    de = pf->GetParam(0);
    a  = pf->GetParam(1);
    r0 = pf->GetParam(2);

    double e1 = exp((1.0-rr/r0)*a/2.0);
    double e2 = e1*e1;

    val  = de*(e2-2.0*e1);
    grad = -de*a/r0*(e2-e1);
  }


// Poly_5 ------------------------------------------------------
  static void Poly_5(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double  r1 = rr - 1.0;
    double dr1 = r1*r1;
    val =     pf->GetParam(0)
        + 0.5*pf->GetParam(1)*dr1
        +     pf->GetParam(2)*r1*dr1
        +     pf->GetParam(3)*dr1*dr1
        +     pf->GetParam(4)*r1*dr1*dr1;
    grad =    pf->GetParam(1)*r1
        + 3.0*pf->GetParam(2)*dr1
        + 4.0*pf->GetParam(3)*r1*dr1
        + 5.0*pf->GetParam(4)*dr1*dr1;
  }

}; // End of Class ===========================================//


