#include <iostream>
#include <math.h>
#include <stdio.h>
#include <tcl.h>
#include <tk.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <blt.h>
#include <time.h>
#include "sample.h"
#include "bilayer.h"
#include "g2model.h"
#define PI 3.1415926535897932
#define Tnum 70
#define Pnum 18
#define SNUM 51

using std::cout; using std::endl;

extern "C" int Blt_Init(Tcl_Interp *interp);
inline double sinc(double a){return(sin(a)/a);}
inline double safesinc(double a){return fabs(a)>1e-10?sin(a)/a:1;}
inline double myrand(){return(rand()%2==1?-rand()/32768.:rand()/32768.);}

typedef double (*Func)(double);

int program_mode=3;
int scatt_mode=0;
// numerical
// 1: the program computes F(q) by numerical Fourier transform over rho(z). 
//    , which is very slow
int numerical=0; 
int amplitudes();
// normal_mode
// 0: pre 2010 with errors being scaled* 
// 1: post 2010 with errors being absolute
// 2: same as 0 except cons = 1 instead of the weird (non)normalization
int normal_mode = 2; 
sample sp[SNUM];
G2model g2;
Tcl_Interp *interp;
char cmd[1024];
const char *colors[] = {"#000000", "#ff0000", "#00ff00", "#0000ff", "#dddd00", 
                  "#00ffff", "#ff00ff", "#999999", "#ff0099", "#9900ff",
                  "#0099ff", "#ff0000", "#00ff00", "#0000ff", "#ffff00", 
                  "#00ffff",
		              "#ff00ff","#999999","#ff0099","#9900ff","#0099ff",
		              "#ff0000","#00ff00","#0000ff","#ffff00","#00ffff",
		              "#ff00ff","#999999","#ff0099","#9900ff","#0099ff",
		              "#ff0000","#00ff00","#0000ff","#ffff00","#00ffff",
		              "#ff00ff","#999999","#ff0099","#9900ff","#0099ff",
		              "#ff0000","#00ff00","#0000ff","#ffff00","#00ffff",
		              "#ff00ff","#999999","#ff0099","#9900ff","#0099ff"};
const char *colorp[] = {"#ffffff", "#ffff00"};
const char *whiteclr = "#ffffff";

Blt_Vector *xmr, *ymCG,*ymPh,*ymCh, *ymc1, *ymc2, *ymc3,*ymw, *ympep;
Blt_Vector *xmr_EDP_v, *ymCG_EDP_v,*ymPh_EDP_v,*ymCh_EDP_v, *ymc1_EDP_v, *ymc2_EDP_v, *ymc3_EDP_v,*ymw_EDP_v,*ympep_EDP_v;
double xmr_EDP[2001],ymCG_EDP[2001],ymPh_EDP[2001],ymCh_EDP[2001],ymc1_EDP[2001],ymc2_EDP[2001],ymc3_EDP[2001],ymw_EDP[2001], ympep_EDP[2001];

long idum;
const char *ctfgraph[2] = {".t.cN", ".t.cX"};
const char *rhograph[2] = {".t.rN", ".t.rX"};
int activesp=0,constraint[2],loops=10,method=1,as_is=1,direct_err=0;
int xpin[Tnum], xindex[Tnum],rdim;

// Arrays for bound-constrained minimization. Linked to tcl arrays with 
// the same names.
bool hasLowerBound[Tnum];
bool hasUpperBound[Tnum];
double upperBounds[Tnum];
double lowerBounds[Tnum];

int anneal_iter=1000,usesign=/*0;NK=*/1,showsign=/*1;NK=*/0,aline_order=1,aline_on_model=0,redraw=1,init_normalize=0;
int max_peak=2;
int noscale=0;//1=don't scale; 0=scale
double anneal_temp=10,stepsize_integral=0.05;
double range=1,tolerance=1e-8;
double sx[Tnum],*x;
double fract[3];
double chifactorN=1e5;
int NMAX=5000;


/* link tcl variables with C variables */
void linkvar(){
  char varname[64];
  int i;
  sprintf(varname ,"set Tnum %d",Tnum);
  Tcl_EvalEx(interp,varname,-1,TCL_EVAL_GLOBAL);
  sprintf(varname ,"set SNUM %d",SNUM);
  Tcl_EvalEx(interp,varname,-1,TCL_EVAL_GLOBAL);
  
  for (int i = 0; i < Tnum; i++) {
    sprintf(varname, "upperBounds(%d)", i);
    Tcl_LinkVar(interp, varname, (char *)&(upperBounds[i]), TCL_LINK_DOUBLE);
  }
  for (int i = 0; i < Tnum; i++) {
    sprintf(varname, "lowerBounds(%d)", i);
    Tcl_LinkVar(interp, varname, (char *)&(lowerBounds[i]), TCL_LINK_DOUBLE);
  }
  for (int i = 0; i < Tnum; i++) {
    sprintf(varname, "hasLowerBound(%d)", i);
    Tcl_LinkVar(interp, varname, (char *)&(hasLowerBound[i]), TCL_LINK_BOOLEAN);
  }
  for (int i = 0; i < Tnum; i++) {
    sprintf(varname, "hasUpperBound(%d)", i);
    Tcl_LinkVar(interp, varname, (char *)&(hasUpperBound[i]), TCL_LINK_BOOLEAN);
  }
  for(i=0;i<Tnum;i++){
    sprintf(varname,"xpin(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(xpin[i]), TCL_LINK_INT);
  }
  for(i=0;i<SNUM;i++){
    sprintf(varname,"status(%d)",i);
	Tcl_LinkVar(interp,varname,(char *)&(sp[i].frmelm.status), TCL_LINK_INT);
  }
  for(i=0;i<2;i++){
    sprintf(varname,"cnst(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(constraint[i]), TCL_LINK_INT);
  }
  
  Tcl_LinkVar(interp,"init_normalize",(char *)&init_normalize, TCL_LINK_INT);
  Tcl_LinkVar(interp,"loops",(char *)&loops, TCL_LINK_INT);
  Tcl_LinkVar(interp,"range",(char *)&range, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"as_is",(char *)&as_is, TCL_LINK_INT);
  Tcl_LinkVar(interp,"direct_err",(char *)&direct_err, TCL_LINK_INT);
  Tcl_LinkVar(interp,"anneal_temp",(char *)&anneal_temp, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"anneal_iter",(char *)&anneal_iter, TCL_LINK_INT);
  Tcl_LinkVar(interp,"method",(char *)&method, TCL_LINK_INT);
  Tcl_LinkVar(interp,"tolerance",(char *)&tolerance, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"usesign",(char *)&usesign, TCL_LINK_INT);
  Tcl_LinkVar(interp,"showsign",(char *)&showsign, TCL_LINK_INT);
  Tcl_LinkVar(interp,"aline_order",(char *)&aline_order, TCL_LINK_INT);
  Tcl_LinkVar(interp,"aline_on_model",(char *)&aline_on_model, TCL_LINK_INT);
  Tcl_LinkVar(interp,"redraw",(char *)&redraw, TCL_LINK_INT);
  Tcl_LinkVar(interp,"stepsize_integral",(char *)&stepsize_integral, TCL_LINK_DOUBLE);
  for(i=0;i<3;i++){
    sprintf(varname,"fract(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(fract[i]), TCL_LINK_DOUBLE);
  }
  Tcl_LinkVar(interp,"program_mode",(char *)&program_mode, TCL_LINK_INT);
  Tcl_LinkVar(interp,"scatt_mode",(char *)&scatt_mode, TCL_LINK_INT);  
  Tcl_LinkVar(interp,"max_peak",(char *)&max_peak, TCL_LINK_INT);  
  Tcl_LinkVar(interp,"noscale",(char *)&noscale, TCL_LINK_INT);  
  Tcl_LinkVar(interp,"numerical",(char *)&numerical, TCL_LINK_INT);  
  Tcl_LinkVar(interp,"chifactorN",(char *)&chifactorN, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"NMAX",(char *)&NMAX, TCL_LINK_INT);  
  Tcl_LinkVar(interp,"normal_mode",(char *)&normal_mode, TCL_LINK_INT);
  char tmp[20];
  strcpy(tmp, "xmr_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &xmr_EDP_v);
  strcpy(tmp, "ymCG_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ymCG_EDP_v);
  strcpy(tmp, "ymPh_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ymPh_EDP_v);
  strcpy(tmp, "ymCh_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ymCh_EDP_v);
  strcpy(tmp, "ymc1_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ymc1_EDP_v);
  strcpy(tmp, "ymc2_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ymc2_EDP_v);
  strcpy(tmp, "ymc3_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ymc3_EDP_v);
  strcpy(tmp, "ymw_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ymw_EDP_v);
  strcpy(tmp, "ympep_EDP_v");
  Blt_CreateVector(interp, tmp, 2001, &ympep_EDP_v);
  
  Tcl_LinkVar(interp, "carbGlyc(target_c)", (char *)&(g2.carbGlyc.target_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "carbGlyc(target_a)", (char *)&(g2.carbGlyc.target_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "carbGlyc(target_s)", (char *)&(g2.carbGlyc.target_s), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "carbGlyc(tol_c)", (char *)&(g2.carbGlyc.tol_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "carbGlyc(tol_a)", (char *)&(g2.carbGlyc.tol_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "carbGlyc(tol_s)", (char *)&(g2.carbGlyc.tol_s), TCL_LINK_DOUBLE);
  
  Tcl_LinkVar(interp, "phosphate(target_c)", (char *)&(g2.phosphate.target_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "phosphate(target_a)", (char *)&(g2.phosphate.target_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "phosphate(target_s)", (char *)&(g2.phosphate.target_s), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "phosphate(tol_c)", (char *)&(g2.phosphate.tol_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "phosphate(tol_a)", (char *)&(g2.phosphate.tol_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "phosphate(tol_s)", (char *)&(g2.phosphate.tol_s), TCL_LINK_DOUBLE);

  Tcl_LinkVar(interp, "choline(target_c)", (char *)&(g2.choline.target_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "choline(target_a)", (char *)&(g2.choline.target_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "choline(target_s)", (char *)&(g2.choline.target_s), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "choline(tol_c)", (char *)&(g2.choline.tol_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "choline(tol_a)", (char *)&(g2.choline.tol_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "choline(tol_s)", (char *)&(g2.choline.tol_s), TCL_LINK_DOUBLE);
  
  Tcl_LinkVar(interp, "methine(target_c)", (char *)&(g2.methine.target_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methine(target_a)", (char *)&(g2.methine.target_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methine(target_s)", (char *)&(g2.methine.target_s), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methine(tol_c)", (char *)&(g2.methine.tol_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methine(tol_a)", (char *)&(g2.methine.tol_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methine(tol_s)", (char *)&(g2.methine.tol_s), TCL_LINK_DOUBLE);
  
  Tcl_LinkVar(interp, "methyl(target_c)", (char *)&(g2.methyl.target_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methyl(target_a)", (char *)&(g2.methyl.target_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methyl(target_s)", (char *)&(g2.methyl.target_s), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methyl(tol_c)", (char *)&(g2.methyl.tol_c), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methyl(tol_a)", (char *)&(g2.methyl.tol_a), TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "methyl(tol_s)", (char *)&(g2.methyl.tol_s), TCL_LINK_DOUBLE);
}

/* this updates the screen */
void updatelink(){
  int i;
  for(i=0;i<Tnum;i++){
    sprintf(cmd,"x(%d)",i);
    Tcl_UpdateLinkedVar(interp, cmd);
  }
  for(i=0;i<3;i++){
    sprintf(cmd,"fract(%d)",i);
    Tcl_UpdateLinkedVar(interp,cmd);
  }
}

void polint(const double *xa, const double *ya,int n,double x,double *y,double *dy)
{
  /* interpolation from Numerical Recipe */
  int i,m,ns=0;
  double den,dif,dift,ho,hp,w;
  double *c=(double *)malloc(n*sizeof(double)),*d=(double *)malloc(n*sizeof(double));
  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns];
  for (m=1;m<n;m++) {
    for (i=0;i<n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      //cout<<den<<endl;
      if ( (den=ho-hp) == 0.0){
	printf("Error in routine POLINT");
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns<(n-m) ? c[ns] : d[--ns]));
  }
  free(d);
  free(c);
}


double polydispcorr(double q, double Rm, double sigR){
    double zSchulz=Rm*Rm/sigR/sigR-1;
    double sSchulz=(zSchulz+1)/Rm;
    double DISTR=exp(-(zSchulz+1)/2*log(1+4*q*q/sSchulz/sSchulz));
	DISTR*=sSchulz*sSchulz;
	DISTR*=cos((3+zSchulz)*atan(2*q/sSchulz));
	DISTR/=-(4*q*q+sSchulz*sSchulz);
	DISTR+=1;
	DISTR*=(zSchulz+2)*(zSchulz+1)/sSchulz/sSchulz/*/q/q*/;
	DISTR*=78.9568352;//8pi^2
	return DISTR;
}


int NK_polydispout(ClientData clientData, Tcl_Interp *interp,int objc, 
                   Tcl_Obj *const objv[]) 
{
  FILE *fp;
  if(objc < 3) return TCL_OK;
  double Rm, sigR;
  Tcl_GetDoubleFromObj(interp,objv[1], &Rm);
  Tcl_GetDoubleFromObj(interp,objv[2], &sigR);
  fp = fopen("polydisp.dat", "w");
  double step = 1.4/2000;
  for (double q = 0.0001; q < 1.4; q += step){
    fprintf(fp, "%11.4g %11.4g\n", q, polydispcorr(q, Rm, sigR));
  }
  fclose(fp);
  return TCL_OK;
}


int Init_vectors()
{
  /* initialize BLT vectors */
  char vectorname[16];
  sprintf(vectorname,"xmr");
  if (Blt_CreateVector(interp, vectorname, 2048, &xmr) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ymCG");
  if (Blt_CreateVector(interp, vectorname, 2048, &ymCG) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ymPh");
  if (Blt_CreateVector(interp, vectorname, 2048, &ymPh) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ymCh");
  if (Blt_CreateVector(interp, vectorname, 2048, &ymCh) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ymc1");
  if (Blt_CreateVector(interp, vectorname, 2048, &ymc1) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ymc2");
  if (Blt_CreateVector(interp, vectorname, 2048, &ymc2) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ymc3");
  if (Blt_CreateVector(interp, vectorname, 2048, &ymc3) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ymw");
  if (Blt_CreateVector(interp, vectorname, 2048, &ymw) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"ympep");
  if (Blt_CreateVector(interp, vectorname, 2048, &ympep) != TCL_OK) return TCL_ERROR;
  return TCL_OK;
}


double eachmctf(int j, int i, double *chi){
  /* supports calculation of the model FF for output in the report command*/
  double *qs,*fs,*er,scale,temp,tmp,theor;
  if(sp[j].frmelm.status==0) return -1;
  qs=sp[j].Q;
  fs=sp[j].TOT;
  er=sp[j].errf;
  scale=sp[j].SCL;
//  amplitudes();
  theor=fabs(numerical==0?sp[j].q2f(qs[i]):sp[j].q2fnumer(qs[i]));
  tmp= (fs[i]==fabs(fs[i]))?fs[i]:0.0;
  switch (normal_mode) {
	  case 0: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-tmp)) / er[i]; break;
	  case 1: temp = (usesign ? (theor-fs[i]*scale) : (theor-tmp*scale)) / er[i]; break;
	  case 2: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-tmp)) / er[i]; break;
	  default: exit(1);
  }
  *chi=temp*temp*(sp[j].mode=='n'?chifactorN:1);
  return theor;
}


/**************
  Start the fitting   
*****************/

double basechi(int j) {
  /* chisquare of all the data points for all selected
   samples without the penalty terms */
  int i,hn;
  double *qs,*fs,*er,scale,chi=0,rsd=0,temp,tmp,theor;
//  for(j=0;j<SNUM;j++){
    if(sp[j].frmelm.status==0) return 0;
	sp[j].scale2model();
    qs=sp[j].Q;
    fs=sp[j].TOT;
    er=sp[j].errf;
    scale=sp[j].SCL;
    hn=sp[j].hn;
	for(i=1;i<=hn;i++){
		theor=fabs(numerical==0?sp[j].q2f(qs[i]):sp[j].q2fnumer(qs[i]));
	    tmp=(fs[i]==fabs(fs[i]))?fs[i]:0.0;
		switch (normal_mode) {
			case 0: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-tmp)) / er[i]; break;
			case 1: temp = (usesign ? (theor-fs[i]*scale) : (theor-tmp*scale)) / er[i]; break;
			case 2: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-tmp)) / er[i]; break;
			default: exit(1);
		}
		chi+=temp*temp;
		temp = (usesign ? (theor-fs[i]*scale) : (theor-tmp*scale));
        rsd+=temp*temp;
    }
	g2.rsdj=rsd*(sp[j].mode=='n'?chifactorN:1);
	sp[j].basek = chi*(sp[j].mode=='n'?chifactorN:1);
	return sp[j].basek;
}


double basechi_sim(int j) {
  int i,hn;
  double *qs,*fs,*er,scale,chi=0,rsd=0,temp,theor;
//  for(j=0;j<SNUM;j++){
    if(sp[j].frmelm.status==0) return 0;
	sp[j].scale2model();
    qs=sp[j].Q;
    fs=sp[j].TOT;
    er=sp[j].errf;
    scale=sp[j].SCL;
    hn=sp[j].hn;
	for(i=1;i<=hn;i++){
		theor = usesign ? sp[activesp].interpolq2f(qs[i]) : fabs(sp[activesp].interpolq2f(qs[i]));
		switch (normal_mode) {
			case 0: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-fabs(fs[i]))) / er[i]; break;
			case 1: temp = (usesign ? (theor-fs[i]*scale) : (theor-fabs(fs[i])*scale)) / er[i]; break;
			case 2: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-fabs(fs[i]))) / er[i]; break;
			default: exit(1);
		}
		chi+=temp*temp;
		temp = (usesign ? (theor-fs[i]*scale) : (theor-fabs(fs[i])*scale));
        rsd+=temp*temp;
    }
	g2.rsdj=rsd*(sp[j].mode=='n'?chifactorN:1);
	sp[j].basek = chi*(sp[j].mode=='n'?chifactorN:1);
	return sp[j].basek;
}


int amplitudes(){
  fract[0]=fract[2]=0; //hard fix (peptide is in the hydrocarbon region only)
  fract[1]=1;
  g2.DC=g2.parC[0];
  if ( (program_mode==1) || (program_mode==2) ) {
  //use program_mode=2 when not using the peptide gaussian
  //(it will speed up a program little bit)
    g2.A=(g2.VL-g2.VHL)/(g2.DC);
  } else if (program_mode==3) {
    g2.A=(g2.VL+fract[1]*g2.Vpep-g2.VHL)/(g2.DC);
  } else return TCL_ERROR;
  g2.VCG=g2.VHL*g2.RCG;
  g2.VPh=g2.VHL*g2.RPh;
  g2.VCh=g2.VHL-g2.VCG-g2.VPh;
  g2.Vc2=(g2.DC*g2.A-fract[1]*g2.Vpep)/2/(g2.nc2+g2.nc1*g2.parr12+g2.nc3*g2.parr);
  g2.Vc3=g2.parr*g2.Vc2;
  g2.Vc1=g2.parr12*g2.Vc2;

  g2.parCG[1]=(g2.parCG[2]*g2.A==0?0:  g2.VCG/g2.parCG[2]/g2.A);
  g2.parPh[1]=(g2.parPh[2]*g2.A==0?0:  g2.VPh/g2.parPh[2]/g2.A);
  g2.parCh[1]=(g2.parCh[2]*g2.A==0?0:  g2.VCh/g2.parCh[2]/g2.A);
  g2.parc3[1]=(g2.parc3[2]*g2.A==0?0:2.*g2.nc3*g2.Vc3/g2.parc3[2]/g2.A);
  g2.parc1[1]=(g2.parc1[2]*g2.A==0?0:2.*g2.nc1*g2.Vc1/g2.parc1[2]/g2.A);
  g2.parC[1]=/*1*/(g2.DC*g2.A-fract[1]*g2.Vpep)/(g2.DC*g2.A);
  return TCL_OK;
}


/******************************************************************************
  calculate the penalty terms 
******************************************************************************/  
double penalty()
{ 
  double pnty = 0, tmp;
  amplitudes();
  g2.VH=g2.VHL+fract[0]*g2.Vpep;
  g2.dXH=g2.parPh[0]-g2.parCG[0];
  g2.dXH2=g2.parPh[0]-g2.parCh[0];
  g2.dXH3=g2.parCG[0]-g2.DC;
  g2.DHH=2*g2.max();
  g2.DH1=g2.DHH/2-g2.DC;
  g2.DB=2*(g2.VL+g2.Vpep)/g2.A;
  if (g2.twiggl>0) g2.wiggl=g2.wiggeling();
  //g2.DBG=g2.gibbs();

  if(g2.tDC>0){tmp=(g2.sDC-g2.DC)/g2.tDC;pnty+=tmp*tmp;}
  if(g2.tr>0){tmp=(g2.sr-g2.parr)/g2.tr;pnty+=tmp*tmp;}
  if(g2.tr12>0){tmp=(g2.sr12-g2.parr12)/g2.tr12;pnty+=tmp*tmp;}
  if(g2.tdXH>0){tmp=(g2.sdXH-g2.dXH)/g2.tdXH;pnty+=tmp*tmp;}
  if(g2.tdXH2>0){tmp=(g2.sdXH2-g2.dXH2)/g2.tdXH2;pnty+=tmp*tmp;}
  if(g2.tdXH3>0){tmp=(g2.sdXH3-g2.dXH3)/g2.tdXH3;pnty+=tmp*tmp;}
  if(g2.tDH1>0){tmp=(g2.DH1R-g2.DH1)/g2.tDH1;pnty+=tmp*tmp;}
  if(g2.tSIGC>0){tmp=(g2.sSIGC-g2.parC[2])/g2.tSIGC;pnty+=tmp*tmp;}
  if(g2.tRCG>0){tmp=(g2.sRCG-g2.RCG)/g2.tRCG;pnty+=tmp*tmp;}
  if(g2.tRPh>0){tmp=(g2.sRPh-g2.RPh)/g2.tRPh;pnty+=tmp*tmp;}
  if(g2.twiggl>0){tmp=(0.0-g2.wiggl)/g2.twiggl;pnty+=tmp*tmp;}
  //if(g2.tDBG>0){tmp=(g2.DB-g2.DBG)/g2.tDBG;pnty+=tmp*tmp;}

  // Go through Gauss objects' soft constraints.
  // First, update their parameters, which are center, ampl, and sigma.
  g2.carbGlyc.set_params(g2.parCG[0], g2.parCG[1], g2.parCG[2]);
  g2.phosphate.set_params(g2.parPh[0], g2.parPh[1], g2.parPh[2]);
  g2.choline.set_params(g2.parCh[0], g2.parCh[1], g2.parCh[2]);
  g2.methine.set_params(g2.parc1[0], g2.parc1[1], g2.parc1[2]);
  g2.methyl.set_params(g2.parc3[0], g2.parc3[1], g2.parc3[2]);
  // Now, collect the penalty values
  pnty += g2.get_penalty();
  
  return pnty;
}

void updatecolor(const char *widget, int k, const char *color)
{
  //fancy things
  sprintf(cmd, "%s%d configure -bg %s", widget, k, color);
  Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL);
}

void updatecolor2(const char *widget, int k, const char *tail, const char *color)
{
  //fancy things
  sprintf(cmd, "%s%d%s configure -bg %s", widget, k, tail, color);
  Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL);
}


/****************************************************************************** 
return chisquare + penalty terms
This is the function that gets minimized by the amoeba function.

*pars: 1D array that holds the values of free parameters with which the current
       function evaluation is being done.
******************************************************************************/
double vec(double *pars)
{  
  for (unsigned int i = 0; i < rdim; i++) {
    // Transform internal variable to external one, then store in the x array.
    // If bounds are not specified for an input parameter, int_to_ext returns 
    // the input value.
    x[xindex[i]] = int_to_ext(xindex[i], pars[i]);
  }
  
  // Make sure the following parameters take on positive values.
  x[0]=fabs(x[0]);x[2]=fabs(x[2]);
  x[3]=fabs(x[3]);x[5]=fabs(x[5]);
  x[6]=fabs(x[6]);x[8]=fabs(x[8]);
  x[9]=fabs(x[9]);x[11]=fabs(x[11]);
  x[12]=fabs(x[12]);x[14]=fabs(x[14]);
  x[15]=fabs(x[15]);x[16]=fabs(x[16]);

  // What is this?
  if (x[4]*x[5]==0) {x[3]=x[0];}

  // Get the penalty due to soft constraints
  g2.penalty = penalty();
  
  g2.basekN = 0, g2.basekX = 0;
  g2.rsdtN = 0, g2.rsdtX = 0;
  for (unsigned int j = 0; j < SNUM; j++) {
	  if (sp[j].frmelm.status==0) continue;
  	if (sp[j].mode=='n') {
  	  g2.basekN+=basechi(j); 
  	  g2.rsdtN+=g2.rsdj;
  	}
    else if (sp[j].mode=='x') {
      g2.basekX += basechi(j); 
      g2.rsdtX += g2.rsdj;
    }
  }
  g2.basek = g2.basekN + g2.basekX;
  g2.rsd = g2.rsdtN + g2.rsdtX;
  return g2.basek + g2.penalty;
}


int YF_display(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // update screen after fitting
  updatelink();
  int i;
  for(i=0;i<SNUM;i++){
    if(sp[i].frmelm.status) sp[i].redrawfrm();
  }
  Tcl_EvalEx(interp,"mdplot",-1,TCL_EVAL_GLOBAL);
  return TCL_OK;
}
int YF_spinfo(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  //print out some fitting setup/result information
  int j,n=0,h=0,c=0,f=Tnum;
  if(g2.tDC>0) c++;
  if(g2.tr>0) c++;
  if(g2.tr12>0) c++;
  if(g2.tdXH>0) c++;
  if(g2.tdXH2>0) c++;
  if(g2.tdXH3>0) c++;
  if(g2.tDH1>0) c++;
  if(g2.tSIGC>0) c++;
  if(g2.tRCG>0) c++;
  if(g2.tRPh>0) c++;

  for(j=0;j<SNUM;j++)	if(sp[j].frmelm.status) {n++, h+=sp[j].hn;}
  for(j=0;j<Tnum;j++)		if(xpin[j]==1) {f--;}
  sprintf(cmd,"MSD=%5lf X2=%5lf X2red=%5lf PT2=%5lf PX2tot=%5lf N=%d v=%d c=%d", g2.rsd/(h-n-f-1+c), g2.basek, g2.basek/(h-n-f-1+c), g2.penalty, g2.basek+ g2.penalty,n,h,c);
  Tcl_SetResult(interp,cmd,TCL_STATIC);
  return TCL_OK;
}
int YF_spinfoX(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  //print out some fitting setup/result information
  sprintf(cmd,"MSDfractX=%5lf X2fractX=%5lf", g2.rsdtX/g2.rsd, g2.basekX/g2.basek);
  Tcl_SetResult(interp,cmd,TCL_STATIC);
  return TCL_OK;
}
int YF_spinfoN(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  //print out some fitting setup/result information
  sprintf(cmd,"MSDfractN=%5lf X2fractN=%5lf", g2.rsdtN/g2.rsd, g2.basekN/g2.basek);
  Tcl_SetResult(interp,cmd,TCL_STATIC);
  return TCL_OK;
}


/******************************************************************************
start the amoeba fitting driver
Pnum: 18
Tnum: 70
xpin: an array that tells whether a parameter is free or fixed.
xindex: an array for the indices of free parameters
x: an array for parameter values, linked to a tcl array named x. 
   See also globalVariables.tcl for more details.
sx[i]:
ran1:
idum:
******************************************************************************/
int YF_amoeba(ClientData clientData, Tcl_Interp *interp, int objc, 
              Tcl_Obj *const objv[]) 
{
  static double *p[Pnum+1], tmpp[Pnum+1], y[Pnum+1];
  double tsx[Pnum], tta, tmp;
  int i, j;
  static char malloced = 0;
  
  if (!malloced) {
    for (i = 0; i < Pnum + 1; i++) p[i] = (double *)malloc(sizeof(double)*Pnum);
    malloced = 1;
  }

  // Get the indices for free parameters
  for (i = 0, rdim = 0; i < Tnum; i++) {
    if(xpin[i]) continue;
    xindex[rdim++] = i;
  }

  // Initialize 2D array p, which determines the initial parameter values
  for (i = 0; i < rdim+1; i++) {
    for (j = 0; j < rdim; j++) p[i][j] = x[xindex[j]];
  }
  // What does this line do?
  for (i = 0; i < rdim; i++) tsx[i] = sx[xindex[i]];
  
  for(i = 1; i < rdim; i++) {
    tta = ran1(&idum) * 2 * PI;
    tmp = tsx[i-1]*cos(tta) - tsx[i]*sin(tta);
    tsx[i] = tsx[i-1]*sin(tta) + tsx[i]*cos(tta);
    tsx[i-1] = tmp;
  }
  
  for (i = 0; i < rdim; i++) p[i][i] += tsx[i];
  
  for (j = 0; j < rdim; j++) tmpp[j] = x[xindex[j]];
  
  // Transform external variables to internal ones
  for (int i = 0; i < rdim+1; i++) {
    for (int j = 0; j < rdim; j++) p[i][j] = ext_to_int(xindex[j], p[i][j]);
  }
  
  // p: initial values
  // y: probably an array that gets used internally by amoeba
  // rdim: the number of free parameters
  // tolerance:
  // vec: the actual function that gets minimized
  if (amoeba(p, y, rdim, tolerance, vec) >= NMAX) {
    // Looks like amoeba did not find the best fit in enough iteration, so
    // put elements in x back to their original values
    for (j = 0; j < rdim; j++) x[xindex[j]] = tmpp[j];
	  NMAX+=5000;
	  return TCL_OK;
  } else {
    // Looks like p[0][j] contains the best fit values
    for (j = 0; j < rdim; j++) x[xindex[j]] = int_to_ext(xindex[j], p[0][j]);
  }

  x[0]=fabs(x[0]);x[2]=fabs(x[2]);
  x[3]=fabs(x[3]);x[5]=fabs(x[5]);
  x[6]=fabs(x[6]);x[8]=fabs(x[8]);
  x[9]=fabs(x[9]);x[11]=fabs(x[11]);
  x[12]=fabs(x[12]);x[14]=fabs(x[14]);
  x[15]=fabs(x[15]);x[16]=fabs(x[16]);
  if(x[0]>x[3]) {tmp=x[0];x[0]=x[3];x[3]=tmp;}
  if (x[4]*x[5]==0) {x[3]=x[0];}

  int n=0,h=0,c=0,f=Tnum;
  if(g2.tDC>0) c++;
  if(g2.tr>0) c++;
  if(g2.tr12>0) c++;
  if(g2.tdXH>0) c++;
  if(g2.tdXH2>0) c++;
  if(g2.tdXH3>0) c++;
  if(g2.tDH1>0) c++;
  if(g2.tSIGC>0) c++;
  if(g2.tRCG>0) c++;
  if(g2.tRPh>0) c++;

  for(j=0;j<SNUM;j++)	if(sp[j].frmelm.status) {n++, h+=sp[j].hn;}
  for(j=0;j<Tnum;j++)		if(xpin[j]==1) {f--;}
  sprintf(cmd,"set msga \"PX2tot=%g\"",g2.basek+g2.penalty);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgb \"PT2=%g\"",g2.penalty);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgc \"X2red=%g\"",g2.basek/(h-n-f-1+c));
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgd \"MSD=%g\"",g2.rsd/(h-n-f-1+c));
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msge \"X2fractN=%g\"",g2.basekN/g2.basek);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgf \"MSDfractN=%g\"",g2.rsdtN/g2.rsd);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgg \"X2fractX=%g\"",g2.basekX/g2.basek);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgh \"MSDfractX=%g\"",g2.rsdtX/g2.rsd);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL); 
// information about chi-squares is for a finnished iteration
// meaning that it is for a previous set of parameters
  if(redraw){
    sprintf(cmd,"display");
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  }
  Tcl_EvalEx(interp,"update",-1,TCL_EVAL_GLOBAL);
  return TCL_OK;
}

int YF_export(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // support the 'export' command, which write all the *.mrho, *.mctf *.rho ... files
  FILE *fp;
  char *str;
  double *qs, *fs, *er;
  int i, j, nv;
  if(objc<3) return TCL_OK;
  str = Tcl_GetString(objv[2]);
  double mctf, chi, scale;
  for(j=0;j<SNUM;j++) {
    if(sp[j].frmelm.status==0) continue;
    strcpy(cmd,str);
    strcat(cmd,".fmf");
    char istr[2];
    sprintf(istr,"%d",sp[j].number);
    fp=fopen(strcat(cmd,istr),"w");
    qs=(sp[j].frmelm.xv)->valueArr;
    fs=(sp[j].frmelm.yv)->valueArr;
    nv=(sp[j].frmelm.xv)->numValues;
    er=(sp[j].errf);
	  switch (normal_mode) {
		  case 0: scale=(sp[j].SCL); break;
		  case 1: scale=1; break;
		  case 2: scale=(sp[j].SCL); break;
		  default: exit(1);
	  }
	  double *err = new double[nv];
    for(i=1;i<nv;i++){

      mctf=eachmctf(j,i,&chi);
	  err[i]=er[i]*scale/(sp[j].mode=='n'?sqrt(chifactorN):1);

	  fprintf(fp,"%11.4g %11.4g %11.4g %11.4g %11.4g %11.4g\n",qs[i],fs[i],mctf,mctf-fs[i],err[i],chi);
    }
    fclose(fp);

    strcpy(cmd,str);
	  if (sp[j].mode=='x') strcat(cmd,".xff");
	  if (sp[j].mode=='n') strcat(cmd,".nff");
    sprintf(istr,"%d",sp[j].number);
    fp=fopen(strcat(cmd,istr),"w");
	  strcpy(cmd,Tcl_GetVar(interp,"StE_header",1));
	  if (strlen(cmd)!=0) fprintf(fp,"%s\n",cmd);
	  if (sp[j].contrast>=0) {
		  sprintf(cmd,"SIMtoEXP_contrast %c %d %g",sp[j].mode,sp[j].contrast,sp[j].rw);
		  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
		  strcpy(cmd,Tcl_GetVar(interp,"StE_contrast",1));
		  if (strlen(cmd)!=0) fprintf(fp,"%s\n",cmd);
	  }
	  fprintf(fp,"q           |F(q)|      deltaF\n");
    for(i=1;i<nv;i++){
	    fprintf(fp,"%11.4g %11.4g %11.4g\n",qs[i],fs[i],err[i]);
    }
    fclose(fp);

	  delete [] err;

    strcpy(cmd,str);
    strcat(cmd,".mctf");
    sprintf(istr,"%d",sp[j].number);
    fp=fopen(strcat(cmd,istr),"w");
    qs=(sp[j].frmelm.xmv)->valueArr;
    fs=(sp[j].frmelm.ymv)->valueArr;
    nv=(sp[j].frmelm.xmv)->numValues;
    for(i=0;i<nv;i++){
      fprintf(fp,"%11.4g %11.4g\n",qs[i],fs[i]);
    }
    fclose(fp);

	  strcpy(cmd,str);
    strcat(cmd,".sdp");
    sprintf(istr,"%d",sp[j].number);
    fp=fopen(strcat(cmd,istr),"w");
    qs=(sp[j].frmelm.xmsp)->valueArr;
    fs=(sp[j].frmelm.ymsp)->valueArr;
    nv=(sp[j].frmelm.xmsp)->numValues;
	  double *fCG=ymCG->valueArr;
	  double *fPh=ymPh->valueArr;
	  double *fCh=ymCh->valueArr;
	  double *fc2=ymc2->valueArr;
	  double *fc1=ymc1->valueArr;
	  double *fc3=ymc3->valueArr;
	  double *fw=ymw->valueArr;
	  double *fpep=ympep->valueArr;
    for(i=0;i<nv;i++){
		  fprintf(fp,"%11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g\n",qs[i],sp[j].eCG/g2.VCG*fCG[i],sp[j].ePh/g2.VPh*fPh[i],sp[j].eCh/g2.VCh*fCh[i],sp[j].ec2/g2.Vc2*fc2[i],sp[j].ec1/g2.Vc1*fc1[i],sp[j].ec3/g2.Vc3*fc3[i],sp[j].rw*fw[i],sp[j].epep/g2.Vpep*fpep[i],fs[i]);
    }
    fclose(fp);
  }
  qs=xmr->valueArr;
  double *fCG=ymCG->valueArr;
  double *fPh=ymPh->valueArr;
  double *fCh=ymCh->valueArr;
  double *fc2=ymc2->valueArr;
  double *fc1=ymc1->valueArr;
  double *fc3=ymc3->valueArr;
  double *fw=ymw->valueArr;
  double *fpep=ympep->valueArr;
  nv=xmr->numValues;
  strcpy(cmd,str);
  fp=fopen(strcat(cmd,".mrho"),"w");
  for(i=0;i<nv;i++){
    fprintf(fp,"%11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g\n",qs[i],/*fs[i],*/fCG[i],fPh[i],fCh[i],fc2[i],fc1[i],fc3[i],fw[i],fpep[i]);
  }
  fclose(fp);
  return TCL_OK;
}

int MR_calculateEDP(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[])
{
  printf("Hello\n");
  // Attain probabilities
  //double *qs=xmr->valueArr;
  double *pCG=ymCG->valueArr;
  double *pPh=ymPh->valueArr;
  double *pCh=ymCh->valueArr;
  double *pc2=ymc2->valueArr;
  double *pc1=ymc1->valueArr;
  double *pc3=ymc3->valueArr;
  double *pw=ymw->valueArr;
  double *ppep=ympep->valueArr;
  int nv=xmr->numValues;
  int i=0;int j=0;
  
  // Find current sample
  for(j=0;j<SNUM;j++){
    if(sp[j].frmelm.status==0) {continue;}
    else {break;}
   }
  // Calculate EDP
   for(i=0;i<nv;i++)
  	{       
		ymCG_EDP[i]=(sp[j].eCG/g2.VCG)*(pCG[i]);
		ymPh_EDP[i]=(sp[j].ePh/g2.VPh)*(pPh[i]);
		ymCh_EDP[i]=(sp[j].eCh/g2.VCh)*(pCh[i]);

		ymc2_EDP[i]=(sp[j].ec2/g2.Vc2)*(pc2[i]);
		ymc1_EDP[i]=(sp[j].ec1/g2.Vc1)*(pc1[i]);
		ymc3_EDP[i]=(sp[j].ec3/g2.Vc3)*(pc3[i]);

		ympep_EDP[i]=(sp[j].epep/g2.Vpep)*(ppep[i]);
		ymw_EDP[i]=pw[i]*(0.333);
	}
    // Transfer to Blt vector format
    Blt_ResetVector(ymCG_EDP_v,ymCG_EDP, nv, nv, TCL_VOLATILE);
    Blt_ResetVector(ymPh_EDP_v,ymPh_EDP, nv, nv, TCL_VOLATILE);
    Blt_ResetVector(ymCh_EDP_v,ymCh_EDP, nv, nv, TCL_VOLATILE);
    Blt_ResetVector(ymc2_EDP_v,ymc2_EDP, nv, nv, TCL_VOLATILE);
    Blt_ResetVector(ymc1_EDP_v,ymc1_EDP, nv, nv, TCL_VOLATILE);
    Blt_ResetVector(ymc3_EDP_v,ymc3_EDP, nv, nv, TCL_VOLATILE);
    Blt_ResetVector(ymw_EDP_v,ymw_EDP, nv, nv, TCL_VOLATILE);
    Blt_ResetVector(ympep_EDP_v,ympep_EDP, nv, nv, TCL_VOLATILE);

   return TCL_OK;

}


int YF_deletesamples(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  int i;
  for(i=0;i<SNUM;i++){
    sp[i].frmelm.status=0;
  }
  return TCL_OK;
}


int YF_colorp(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  //fancy things
  double *f=sp[activesp].F;
  int hnt = sp[activesp].hn;
  const char *tmp; 
  for(int i = 1; i <= hnt; i++) {
    tmp = (f[i] > 0)? colorp[1]: colorp[0];
    updatecolor(".f.p.", i, tmp);
  }
  return TCL_OK;
}


int YF_parameter(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // handles the 'parameter' command in smp files
  int num;
  Tcl_GetIntFromObj(interp,objv[1],&num);
  sp[num].assignvalues(objc-2,objv+2);
  return TCL_OK;
}


int NK_setsigns(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // triggers the sign assignment command
  int num;
  Tcl_GetIntFromObj(interp,objv[1],&num);
  sp[num].assignsigns(objc-2,objv+2);
  sp[num].redrawfrm();
  return TCL_OK;
}


int YF_newsample(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // initialize a new sample structure, partially implemented 'samplist' command
  int num;
  Tcl_GetIntFromObj(interp,objv[1],&num);
  sp[num].number=num;
  sp[num].initialize(Tcl_GetString(objv[2]));
  return TCL_OK;
}


int YF_plot(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // plot the sample datasets
  int num;
  Tcl_GetIntFromObj(interp, objv[2], &num);
  if (strcmp(Tcl_GetString(objv[1]), "f") != 0) {
    sp[num].drawfrm();
	  sp[num].SCL=1.;
	  const char *tmp = (sp[num].frmelm.status == 0)? whiteclr: colors[num];
    updatecolor2(".f.t.t.", num, ".2", tmp);
  }
  return TCL_OK;
}


int YF_modelplot(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // plot the model
  double *qs,*fs,*fCG,*fPh,*fCh,*fc2,*fc1,*fc3,*fw,*fpep;
  int points;
//  amplitudes();
  for(int j=0;j<SNUM;j++){
    if((sp[j].frmelm.status==0)/*||j==activesp*/) continue;
    qs=sp[j].frmelm.xmv->valueArr;
    fs=sp[j].frmelm.ymv->valueArr;
    points=sp[j].cft(qs,fs);
	Blt_ResetVector(sp[j].frmelm.xmv,qs, points, points, TCL_VOLATILE);
	Blt_ResetVector(sp[j].frmelm.ymv,fs, points, points, TCL_VOLATILE);

	qs=sp[j].frmelm.xmsp->valueArr;
    fs=sp[j].frmelm.ymsp->valueArr;
    points=sp[j].sdp(qs,fs);
	Blt_ResetVector(sp[j].frmelm.xmsp,qs, points, points, TCL_VOLATILE);
	Blt_ResetVector(sp[j].frmelm.ymsp,fs, points, points, TCL_VOLATILE);

    sprintf(cmd,"raiseq %s %d",sp[j].frmelm.graph,sp[j].frmelm.id);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);    
  }
  qs=xmr->valueArr;
//  fs=ymdr->valueArr;
  fCG=ymCG->valueArr;
  fPh=ymPh->valueArr;
  fCh=ymCh->valueArr;
  fc2=ymc2->valueArr;
  fc1=ymc1->valueArr;
  fc3=ymc3->valueArr;
  fw=ymw->valueArr;
  fpep=ympep->valueArr;
  points=g2.rho(qs,/*fs,*/fCG,fPh,fCh,fc2,fc1,fc3,fw,fpep);
  Blt_ResetVector(xmr,qs, points, 2048, TCL_VOLATILE);
//  Blt_ResetVector(ymdr,fs, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ymCG,fCG, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ymPh,fPh, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ymCh,fCh, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ymc2,fc2, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ymc1,fc1, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ymc3,fc3, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ymw,fw, points, 2048, TCL_VOLATILE);
  Blt_ResetVector(ympep,fpep, points, 2048, TCL_VOLATILE);
  return TCL_OK;
}


int NK_calcmodel(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  //calculates the model and scales the sample to model
  int i,j;
  double tmp;
  //if (g2.twiggl<=0 && g2.tDBG>0) g2.tDBG*=(-1.);
  g2.parCG[0]=fabs(g2.parCG[0]); g2.parCG[2]=fabs(g2.parCG[2]);
  g2.parPh[0]=fabs(g2.parPh[0]); g2.parPh[2]=fabs(g2.parPh[2]);
  g2.parCh[0]=fabs(g2.parCh[0]); g2.parCh[2]=fabs(g2.parCh[2]);
  g2.parC[0]=fabs(g2.parC[0]); g2.parC[2]=fabs(g2.parC[2]);
  g2.parc1[0]=fabs(g2.parc1[0]); g2.parc1[2]=fabs(g2.parc1[2]);
  g2.parc3[0]=fabs(g2.parc3[0]); g2.parc3[2]=fabs(g2.parc3[2]);
  if (g2.parCG[0]>g2.parPh[0]) {tmp=g2.parCG[0];g2.parCG[0]=g2.parPh[0];g2.parPh[0]=tmp;}
  if (g2.parPh[1]*g2.parPh[2]==0) {g2.parPh[0]=g2.parCG[0];}
  int n=0,h=0,c=0,f=Tnum;
  if(g2.tDC>0) c++;
  if(g2.tr>0) c++;
  if(g2.tr12>0) c++;
  if(g2.tdXH>0) c++;
  if(g2.tdXH2>0) c++;
  if(g2.tdXH3>0) c++;
  if(g2.tDH1>0) c++;
  if(g2.tSIGC>0) c++;
  if(g2.tRCG>0) c++;
  if(g2.tRPh>0) c++;
  //if(g2.tDBG>0) c++;
  for(j=0;j<SNUM;j++)	if(sp[j].frmelm.status) {n++, h+=sp[j].hn;}
  for(i=0;i<Tnum;i++)		if(xpin[i]==1) {f--;}
  
  g2.penalty=penalty();
  g2.basekN=0;
  g2.basekX=0;
  g2.rsdtN=0;
  g2.rsdtX=0;
  for (j=0;j<SNUM;j++){
	  if (sp[j].frmelm.status==0) continue;
	  if (sp[j].mode=='n') {g2.basekN+=basechi(j); g2.rsdtN+=g2.rsdj;}
	  else if (sp[j].mode=='x') {g2.basekX+=basechi(j); g2.rsdtX+=g2.rsdj;}
  }
  g2.basek=g2.basekN+g2.basekX;
  g2.rsd=g2.rsdtN+g2.rsdtX;
  sprintf(cmd,"set msga \"PX2tot=%g\"",g2.basek+g2.penalty);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgb \"PT2=%g\"",g2.penalty);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgc \"X2red=%g\"",g2.basek/(h-n-f-1+c));
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  //msgi added by KA
  sprintf(cmd,"set msgi \"X2=%g\"",g2.basek);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgd \"MSD=%g\"",g2.rsd/(h-n-f-1+c));
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msge \"X2fractN=%g\"",g2.basekN/g2.basek);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgf \"MSDfractN=%g\"",g2.rsdtN/g2.rsd);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgg \"X2fractX=%g\"",g2.basekX/g2.basek);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msgh \"MSDfractX=%g\"",g2.rsdtX/g2.rsd);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL); 

  return TCL_OK;
}


int NK_scale2sim(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  int i;
  Tcl_GetIntFromObj(interp,objv[1],&activesp);
  double rsd_simi=0,basechi_simi=0;
  for(i=0;i<SNUM;i++){
    if((sp[i].frmelm.status==0)||i==activesp) continue;
    sp[i].scale2sim();
	basechi_simi+=basechi_sim(i);
	rsd_simi+=g2.rsdj;
  }
  sprintf(cmd,"set msgd \"ctf_rsd=%g\"",rsd_simi);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  sprintf(cmd,"set msg6 \"ctf_chi=%g\"",basechi_simi);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  return TCL_OK;
}
int NK_calcFourier(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  //calculates the discrete Fourier transform
  double *qs,*fs,*zP,*rP;
  int points;
  for(int j=0;j<SNUM;j++){
    if((sp[j].frmelm.status==0)/*||j==activesp*/) continue;
	zP=sp[j].frmelm.xmsp->valueArr;
    rP=sp[j].frmelm.ymsp->valueArr;
	points=sp[j].RHOdft(zP,rP);
	Blt_ResetVector(sp[j].frmelm.xmsp,zP, points, points, TCL_VOLATILE);
	Blt_ResetVector(sp[j].frmelm.ymsp,rP, points, points, TCL_VOLATILE);

	qs=sp[j].frmelm.xmv->valueArr;
    fs=sp[j].frmelm.ymv->valueArr;
    points=sp[j].dft(qs,fs,zP,rP);
	Blt_ResetVector(sp[j].frmelm.xmv,qs, points, points, TCL_VOLATILE);
	Blt_ResetVector(sp[j].frmelm.ymv,fs, points, points, TCL_VOLATILE);

    sprintf(cmd,"raiseq %s %d",sp[j].frmelm.graph,sp[j].frmelm.id);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);    
  }
  return TCL_OK;
}


int YF_getreport(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // supports the 'report' command
  int num,i;
  char *str;
  double mctf,tmp,scale,chi;
  Tcl_GetIntFromObj(interp,objv[1],&num);
  str=Tcl_GetString(objv[2]);
  scale=sp[num].SCL;
  sprintf(cmd,"%s insert end \"X2 : %11.4g / %d = %11.4g\n\"",str,sp[num].basek,sp[num].hn,sp[num].basek/sp[num].hn);
  Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  for(i=1;i<=sp[num].hn;i++){
    //sprintf(cmd,"%s insert end \"%d\t%f\t%f\t%f\n\"",str,i,sp[num].F[i],sp[num].errf[i],eachchi(num,i));
    mctf=eachmctf(num,i,&chi);
    if (usesign)  tmp=sp[num].TOT[i];
     else tmp=(sp[num].TOT[i]==fabs(sp[num].TOT[i]))?sp[num].TOT[i]:0.0;
    double error=sp[num].errf[i]/(sp[num].mode=='n'?sqrt(chifactorN):1);
	switch (normal_mode) {
		case 0: error*=scale; break;
		case 1: error*=1;     break;
		case 2: error*=scale; break;
		default: exit(1);
	}
	sprintf(cmd,"%s insert end \"%11.4g\t%11.4g\t%11.4g\t%11.4g\t%11.4g\t%11.4g\n\"",str,sp[num].Q[i],tmp*scale,mctf,mctf-tmp*scale,error,chi);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
  }
  return TCL_OK;
}


int YF_err_report(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]) {
  // output error information to a file
  FILE *fp;
  int i,j,hn;
  double *qs,*fs,*er,scale,chi=0,temp,tmp,theor;
  if(objc>1) fp=fopen(Tcl_GetString(objv[1]),"w");
  for(j=0;j<SNUM;j++){
    if(sp[j].frmelm.status==0) continue;
    chi=0;
    qs=sp[j].Q;
    fs=sp[j].TOT;
    er=sp[j].errf;
    scale=sp[j].SCL;
    hn=sp[j].hn;
	  for(i=1;i<=hn;i++){
		  theor=fabs(numerical==0?sp[j].q2f(qs[i]):sp[j].q2fnumer(qs[i]));
	    tmp=(fs[i]==fabs(fs[i]))?fs[i]:0.0;
		  switch (normal_mode) {
			  case 0: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-tmp)) / er[i]; break;
			  case 1: temp = (usesign ? (theor-fs[i]*scale) : (theor-tmp*scale)) / er[i]; break;
			  case 2: temp = (usesign ? (theor/scale-fs[i]) : (theor/scale-tmp)) / er[i]; break;
			  default: exit(1);
		  }
		  chi+=temp*temp;
    }
    if(objc>1) fprintf(fp,"%d: %f\n",j,chi);
    else printf("%d: %f\n",j,chi);
  }
  if(objc>1) fclose(fp);
  return TCL_OK;
}


int main(int argc, char *argv[])
{
  interp = Tcl_CreateInterp();
  Tk_MainEx(argc, argv, Tcl_AppInit, interp);
  return(0);
}


int Tcl_AppInit(Tcl_Interp *interp)
{
  /* command registration */
  int i;
  if (Tcl_Init(interp) != TCL_OK) {
    printf("%s\n", interp->result);
    return TCL_ERROR;
  } 
  if (Tk_Init(interp) != TCL_OK) { 
    return TCL_ERROR; 
  }
  if (Blt_Init(interp) !=TCL_OK) {
    return TCL_ERROR;
  }
  
  Tcl_CreateObjCommand(interp,"plot", YF_plot, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"newsample", YF_newsample, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"parameter", YF_parameter, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"setsigns", NK_setsigns, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"colorp", YF_colorp, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"export", YF_export, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"calculateEDP", MR_calculateEDP, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"deletesamples", YF_deletesamples, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"mdplot", YF_modelplot, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"calcmodel", NK_calcmodel, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"calcFourier", NK_calcFourier, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"amoeba", YF_amoeba, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"display", YF_display, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"getreport", YF_getreport, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"scale2sim", NK_scale2sim, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"err_report", YF_err_report, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"spinfo", YF_spinfo, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"spinfoX", YF_spinfoX, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"spinfoN", YF_spinfoN, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"polydispout", NK_polydispout, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "setNePep", setNePep, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  
  for(i=0;i<Tnum;i++) {xpin[i]=1;}
  int j[]={0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 19, 20, 21, 22, 23};
  for(i=0;i<Pnum;i++) {xpin[j[i]]=0;}
  for(i=0;i<Pnum;i++) {sx[j[i]]=0.1;}
  Init_vectors();
  g2.init();
  x = g2.parCG;
  linkvar();
  idum=time(NULL);
  Tcl_EvalEx(interp,"source tcl/bilayer.tcl",-1,TCL_EVAL_GLOBAL);
  Tcl_EvalEx(interp,"wm title . \"Bilayer by Norbert.Kucerka@nrc.gc.ca\"",-1,TCL_EVAL_GLOBAL);
  return TCL_OK;
}


/******************************************
amoeba fitting driver from Numerical Recipe
*******************************************/

#define TINY 1.0e-10
//#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define GET_PSUM for (j=0;j<ndim;j++) { for (i=0,sum=0.0;i<mpts;i++) sum += p[i][j]; psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
double amotry(double **p,double *y,double *psum,int ndim, double (*funk)(double *),int ihi, double fac)
{
  int j;
  double fac1,fac2,ytry,*ptry;

  ptry=(double *)malloc(sizeof(double)*ndim);
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  ytry=(*funk)(ptry);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=0;j<ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  free(ptry);
  return ytry;
}


int amoeba(double **p, double *y,int ndim, double ftol, double (*funk)(double *))
{
  int i,j,ilo,ihi,inhi,mpts=ndim+1,nfunk=0;
  double ytry,ysave,sum,rtol,*psum,swap;

  psum=(double *)malloc(sizeof(double)*ndim);
  GET_PSUM
    for(i=0;i<mpts;i++) y[i]=(*funk)(p[i]);
  for (;;) {
    ilo=0;
    ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
    for (i=0;i<mpts;i++) {
      if (y[i] < y[ilo]) ilo=i;
      if (y[i] > y[ihi]) {
	inhi=ihi;
	ihi=i;
      } else if (y[i] > y[inhi] && i!=ihi) inhi=i;
    }
    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
    if (rtol < ftol){
      SWAP(y[0],y[ilo])
	for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
			       break;
    }
    if (nfunk >= NMAX){
      printf("Too many (%d) iterations in AMOEBA\n",NMAX);
      return nfunk;
    }
    ytry=amotry(p,y,psum,ndim,funk,ihi,-ALPHA);
    nfunk++;
    if (ytry <= y[ilo]){
      ytry=amotry(p,y,psum,ndim,funk,ihi,GAMMA);
      nfunk++;
    } else if (ytry >= y[inhi]) {
      ysave=y[ihi];
      ytry=amotry(p,y,psum,ndim,funk,ihi,BETA);
      nfunk++;
      if (ytry >= ysave) {
	for (i=0;i<mpts;i++) {
	  if (i != ilo) {
	    for (j=0;j<ndim;j++) {
	      psum[j]=0.5*(p[i][j]+p[ilo][j]);
	      p[i][j]=psum[j];
	    }
	    y[i]=(*funk)(psum);
	  }
	}
	nfunk += ndim;
	GET_PSUM
	  }
    }
  }
  free(psum);
  return nfunk;
}

#undef ALPHA
#undef BETA
#undef GAMMA
//#undef NMAX
#undef GET_PSUM

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}


double active_interpol(double qq){
  return sp[activesp].interpolq2f(qq);
}


int setNePep(ClientData clientData, Tcl_Interp *interp, int objc, 
             Tcl_Obj *const objv[])
{
  // If head = 1, peptide is assumed to be in the head region. If head != 1, 
  // peptide is assumed to be in the tail region.
  int head;
  Tcl_GetIntFromObj(interp, objv[1], &head);
  cout << "head = " << head << endl;
  // NePep is the number of electrons in the peptide. This value gets assigned 
  // to either eCh or ec1 of sp. eCh is the number of electrons in Choline. ec1 
  // is that in methine (double bond).
  double NePep;
  Tcl_GetDoubleFromObj(interp, objv[2], &NePep);
  cout << "NePep = " << NePep << endl;

  // Go through all samples and find highlighted ones. Then, modify eCh or ec1 
  // of the highlighted samples.
  for (int j = 0; j < SNUM; j++) {
    if (sp[j].frmelm.status == 0) {
      // If a sample is not highlighted in the GUI, status is equal to zero.
      continue;
    }
    cout << "sp[" << j << "]" << ": status = 1" << endl;
    if (head == 1) {
      sp[j].eCh = NePep;
      sp[j].ec1 = 0;
    } else {
      sp[j].ec1 = NePep;
      sp[j].eCh = 0;
    }
  }
  return TCL_OK;
}

// index n and internal value, val
double int_to_ext(unsigned int n, double val)
{
  // Check whether bounds exist for n
  if (hasLowerBound[n] || hasUpperBound[n]) {
    double a, b;
    if (hasLowerBound[n] && hasUpperBound[n]) {
      a = lowerBounds[n];
      b = upperBounds[n];
      return a + (b-a)/2*(sin(val)+1);
    } else if (hasLowerBound[n]) {
      a = lowerBounds[n];
      return a - 1 + sqrt(val*val+1);
    } else if (hasUpperBound[n]) {
      b = upperBounds[n];
      return b + 1 - sqrt(val*val+1);
    }
  }
  
  return val;
}


// index n and external value, val
double ext_to_int(unsigned int n, double val)
{
  // Check whether a bound exists for n
  if (hasLowerBound[n] || hasUpperBound[n]) {
    double a, b;
    if (hasLowerBound[n] && hasUpperBound[n]) {
      a = lowerBounds[n];
      b = upperBounds[n];
      return asin(2*(val-a)/(b-a) - 1);
    } else if (hasLowerBound[n]) {
      a = lowerBounds[n];
      return sqrt((val-a+1)*(val-a+1) - 1);
    } else if (hasUpperBound[n]) {
      b = upperBounds[n];
      return sqrt((b-val+1)*(b-val+1) - 1);
    }
  }
  
  return val;
}
