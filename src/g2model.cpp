#include <tcl.h>
#include <stdio.h>
#include <tk.h>
#include <math.h>

#include "g2model.h"
#include "sample.h"

#define PI 3.1415926535897932
#define Tnum 70
#define Pnum 18
#define SNUM 51

extern G2model g2;
extern int numerical;
extern int usesign;
extern int noscale;
extern char cmd[];
extern Tcl_Interp *interp;
extern char *ctfgraph[];
extern char *rhograph[];
extern int showsign;
extern char *colors[];
extern int direct_err;
extern sample *sp;
extern int max_peak;

void G2model::init(){
  /* initialize parameters for the model */
  nc3=1.;
  linkvar();
}


void G2model::linkvar(){
  /* link tcl variables with C variables */
  char varname[64];
  int i;
  for(i=0;i<Tnum;i++){
    sprintf(varname,"x(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(parCG[i]), TCL_LINK_DOUBLE);
  }
}


double G2model::z2rhoCG(double z){
  /* given z return \rho of the first (CG) gaussian*/
  double t,f,tm;
  t=(parCG[2]==0?0:(z-parCG[0])/parCG[2]);tm=(parCG[2]==0?0:(-z-parCG[0])/parCG[2]);
  f=parCG[1]*(exp(-t*t/2.0)+exp(-tm*tm/2.0))/2.5066283;
  return f;
}


double G2model::z2rhoPh(double z){
  /* given z return \rho of the second (Ph) gaussian*/
  double t,f,tm;
  t=(parPh[2]==0?0:(z-parPh[0])/parPh[2]);tm=(parPh[2]==0?0:(-z-parPh[0])/parPh[2]);
  f=parPh[1]*(exp(-t*t/2.0)+exp(-tm*tm/2.0))/2.5066283;
  return f;
}


double G2model::z2rhoCh(double z){
  /* given z return \rho of the second (Ch) gaussian*/
  double t,f,tm;
  t=(parCh[2]==0?0:(z-parCh[0])/parCh[2]);tm=(parCh[2]==0?0:(-z-parCh[0])/parCh[2]);
  f=parCh[1]*(exp(-t*t/2.0)+exp(-tm*tm/2.0))/2.5066283;
  return f;
}


double G2model::z2rhoC(double z){
  /* given z return \rho of the hydrocarbon region (C) error function*/
  double t,f,tm;
  t=(parC[2]==0?0:(z+parC[0])/parC[2]);tm=(parC[2]==0?0:(z-parC[0])/parC[2]);
  f=1./2.*(erff(t/1.414213562)-erff(tm/1.414213562));	//full calculation
  f*=parC[1];
//  f+=parc[1]/2*(1-erff(tm/1.414213562));
  return f;
}


double G2model::z2rhoc3(double z){
  /* given z return \rho of the third (M) gaussian*/
  double t,f,tm;
  t=(parc3[2]==0?0:(z-parc3[0])/parc3[2]);tm=(parc3[2]==0?0:(-z-parc3[0])/parc3[2]);
  f=parc3[1]*(exp(-t*t/2.0)+exp(-tm*tm/2.0))/2.5066283;
  return f;
}


double G2model::z2rhoc1(double z){
  /* given z return \rho of the third (CH) gaussian*/
  double t,f,tm;
  t=(parc1[2]==0?0:(z-parc1[0])/parc1[2]);tm=(parc1[2]==0?0:(-z-parc1[0])/parc1[2]);
  f=parc1[1]*(exp(-t*t/2.0)+exp(-tm*tm/2.0))/2.5066283;
  return f;
}


double G2model::z2rhoc2(double z){
  /* given z return \rho of the CH2 distribution*/
  double f;
  f=z2rhoC(z)-parC[1]*z2rhoc3(z)-parC[1]*z2rhoc1(z);
  return f;
}


double G2model::z2rhopep(double z){
  /* given z return \rho of the peptide distribution (Pp(z)=cp/cc2*Pc2(z))*/
  double f;
  f=(1-parC[1])/parC[1]*z2rhoc2(z);
  return f;
}


double G2model::z2rhoW(double z){
  /* given z return \rho of the water distribution*/
  double f=1.0;
  rhoc2=z2rhoc2(z);
  f-=z2rhoCG(z);
  f-=z2rhoPh(z);
  f-=z2rhoCh(z);
  f-=1/parC[1]*z2rhoC(z);
  return f;
}


int G2model::rho(double *zP,/*double *rP,*/double *rPCG,double *rPPh,double *rPCh,double *rPc2,double *rPc1,double *rPc3,double *rPW,double *rPpep){
  /* use z2rho to construct the entire EDP of the model */
  int i,PT=1000;
  double STEP,z;
  STEP=D/PT/2.0;
  for(i=0,z=0;i<=PT;i++,z+=STEP){
//    rP[i+PT]=z2rho(z);
    rPCG[i+PT]=z2rhoCG(z);
    rPPh[i+PT]=z2rhoPh(z);
	rPCh[i+PT]=z2rhoCh(z);
    rPc2[i+PT]=z2rhoc2(z);
	rPc1[i+PT]=z2rhoc1(z);
	rPc3[i+PT]=z2rhoc3(z);
    rPW[i+PT]=z2rhoW(z);
	rPpep[i+PT]=z2rhopep(z);
    zP[i+PT]=z;
  }
  for(i=0;i<PT;i++) {
    zP[i]=-zP[2*PT-i];
//    rP[i]=rP[2*PT-i];
    rPCG[i]=rPCG[2*PT-i];
    rPPh[i]=rPPh[2*PT-i];
	rPCh[i]=rPCh[2*PT-i];
    rPc2[i]=rPc2[2*PT-i];
	rPc1[i]=rPc1[2*PT-i];
	rPc3[i]=rPc3[2*PT-i];
    rPW[i]=rPW[2*PT-i];
	rPpep[i]=rPpep[2*PT-i];
  }
  return 2001;
}


double G2model::max(){
  /* find the z position of the maxima of the model */
  /* it is peak-peak distance if max_peak=0
      or HG1-HG1 distance if max_peak=1
	  or HG2-HG2 distance if max_peak=2
	  or HGp-HGp distance if max_peak=3 */
//  if (max_peak==0) max_peak=2;
  switch(max_peak) {
    case 0: {
      for (int j=0;j<SNUM;j++) {
        if (sp[j].frmelm.status==0) continue;
        double step=/*0.01*/twiggl,z=parPh[0],v=sp[j].z2sdp(z),limit=D/2,tmp;
        if(sp[j].z2sdp(z+step)<v) step=-step;
        while(z<limit && z>0){
          z+=step;
          if((tmp=sp[j].z2sdp(z))<v) return z-step;
          else v=tmp;
        }
        return z;
	    }
	    break;
    }
    case 1:
      return parCG[0];
	    break;
    case 2:
      return parPh[0];
	    break;
    case 3:
      return parCh[0];
	    break;
  }
}


double G2model::wiggeling()
{
  double tmp = 0, temp;
  int i;
  double STEP = twiggl, z;
//  STEP=D/PT/2.0;
// speed up the calculation
  if (STEP <= 0) return tmp;
  for(i = 0, z =0; z <= D/2.0 ; i++ , z+=STEP) {
    temp=z2rhoW(z);
    if (temp<0) tmp-=temp;
	  if (rhoc2<0) tmp-=rhoc2;
  }
  return tmp;
}


/******************************************************************************
Collect the penalty term from each Gauss object, sum them up, and return the 
sum.
******************************************************************************/
double G2model::get_penalty()
{
  double tmp, ret = 0;
  
  ret += carbGlyc.penalty();
  ret += phosphate.penalty();
  ret += choline.penalty();
  ret += methine.penalty();
  ret += methyl.penalty();
  
  return ret;    
}


/******************************************************************************

******************************************************************************/
double Gauss::penalty()
{
  double ret = 0, tmp;
  if (tol_c > 0) {
    tmp = (center-target_c) / tol_c;
    ret += tmp*tmp;
  }
  if (tol_a > 0) {
    tmp = (ampl-target_a) / tol_a;
    ret += tmp*tmp;
  }
  if (tol_s > 0) {
    tmp = (sigma-target_s) / tol_s;
    ret += tmp*tmp;
  }
  
  return ret;
}
