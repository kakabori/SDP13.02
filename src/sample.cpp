#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <tk.h>
#include <blt.h>

#include "sample.h"
#include "g2model.h"
#include "bilayer.h"

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
extern char* ctfgraph[];
extern char* rhograph[];
extern int showsign;
extern char *colors[];
extern int direct_err;
extern int normal_mode;
extern double stepsize_integral;

using std::cout; using std::endl;

sample::sample(){
  /* initialization and memory allocation */
  D = 60; P = 1; A = 48; SCL = 1; hn = 2;
  F=(double*)calloc(hn+1, sizeof(double));
  Q=(double*)calloc(hn+1, sizeof(double));
  errf=(double*)calloc(hn+1, sizeof(double));
  TOT=(double*)calloc(hn+1, sizeof(double));
}

sample::~sample(){
  /* memory deallocation */
  free(F);
  free(Q);
  free(errf);
  free(TOT);
}
void sample::scale2model() {
//   obtain a scaling factor to compare the sample with the model 
// by scaling sample to the model (NK in Jan 2010)
  int h;
  double a=0,b=0,c=0,mf,temp;
  if(frmelm.status==0) return;
  for(h=1;h<=hn;h++){    
	if (polydisp==1 && g2.sigR>0 && g2.Rm>50) TOT[h]=F[h]/sqrt(polydispcorr(Q[h],g2.Rm,g2.sigR));
	else TOT[h]=F[h];
	mf=fabs((numerical==0?q2f(Q[h]):q2fnumer(Q[h]))/errf[h]);
    temp= (TOT[h]==fabs(TOT[h])) ? (TOT[h]/errf[h]) : 0.0;
    a+=usesign?(mf*TOT[h]/errf[h]):(mf*temp);
	b+=mf*mf;
	c+=usesign?(TOT[h]/errf[h]*TOT[h]/errf[h]):temp*temp;
  }
  switch (normal_mode) {
	  case 0: SCL = fabs(noscale==1 ? 1 : b/a); break; //cout << "SCL = " << SCL << endl; break;
	  case 1: SCL = fabs(noscale==1 ? 1 : a/c); break; //cout << "SCL = " << SCL << endl; break;
	  case 2: SCL = fabs(noscale==1 ? 1 : b/a); break;
  }
}
void sample::scale2sim() {
  /* obtain a scaling factor to compare the 
	 sample with the simulation */
  int h;
  double a=0,b=0,c=0,mf,temp;
  if(frmelm.status==0) return;
  for(h=1;h<=hn;h++){    
    mf=fabs(active_interpol(Q[h])/errf[h]);
    temp= (TOT[h]==fabs(TOT[h])) ? (TOT[h]/errf[h]) : 0.0;
    a+=usesign?(mf*TOT[h]/errf[h]):(mf*temp);
    b+=mf*mf;
	c+=usesign?(TOT[h]/errf[h]*TOT[h]/errf[h]):temp*temp;
  }
  switch (normal_mode) {
	  case 0: SCL = fabs(noscale==1 ? 1 : b/a); break;
	  case 1: SCL = fabs(noscale==1 ? 1 : a/c); break;
	  case 2: SCL = fabs(noscale==1 ? 1 : b/a); break;
  }
  printf("sample # %d SCL = %lf\n",number,SCL);
}

void sample::redrawfrm() {
  /* update form factors on the screen */
  double *qs,*fs;
  int i;
  if(frmelm.status==0) return;
  qs=frmelm.xv->valueArr;
  fs=frmelm.yv->valueArr;
  if (usesign) {
    for(i=1;i<=hn;i++){
		fs[i]=(TOT[i]*SCL);
	}
  } else {
	  for(i=1;i<=hn;i++){
		  fs[i]=(TOT[i]==fabs(TOT[i]))?((TOT[i])*SCL):0.0 ;
	  }
  }
  if(hn>60){
    Blt_ResetVector(frmelm.xv,qs, hn+1, hn+1, TCL_VOLATILE);
    Blt_ResetVector(frmelm.yv,fs, hn+1, hn+1, TCL_VOLATILE);
  }else{
    Blt_ResetVector(frmelm.xv,qs, hn+1, 32, TCL_STATIC);
    Blt_ResetVector(frmelm.yv,fs, hn+1, 32, TCL_STATIC);
  }
}

void sample::initialize(char *initname){
  /* more initialization of the sample; overwritten by .smp input */
  type=1;hn=0;D=60;P=1;A=48;SCL=1;
  contrast=-1;
  strcpy(name,initname);
  strcpy(color,"#ff0000");
  frmelm.status=0;
  frmelm.id=number;
  sprintf(cmd,"xfrm%d",number);
  Blt_CreateVector(interp, cmd, 32, &(frmelm.xv));
  sprintf(cmd,"yfrm%d",number);
  Blt_CreateVector(interp, cmd, 32, &(frmelm.yv));

  sprintf(cmd,"xmdc%d",number);
  Blt_CreateVector(interp, cmd, 1999, &(frmelm.xmv));
  sprintf(cmd,"ymdc%d",number);
  Blt_CreateVector(interp, cmd, 1999, &(frmelm.ymv));

  sprintf(cmd,"xmsp%d",number);
  Blt_CreateVector(interp, cmd, 2001, &(frmelm.xmsp));
  sprintf(cmd,"ymsp%d",number);
  Blt_CreateVector(interp, cmd, 2001, &(frmelm.ymsp));
}

int sample::assignvalues(int objc, Tcl_Obj *const objv[]){
  /* this processes the 'parameter' command in the smp files 
	That command first reads in the real parameters from one
	'line' (lambda, etc.).  It then keeps going to read in all
	the form factors which are technically in the same computer
	"line"  */
  int i,j;
  char *file;
  file=Tcl_GetString(objv[0]);
  if(beamfile!=NULL) free(beamfile);
  beamfile=(char *)malloc(strlen(file)+1);
  strcpy(beamfile,file);
  if(!strcmp("curved",beamfile)) type=2;
//  if(!strcmp("simul",beamfile)) type=3;
  Tcl_GetDoubleFromObj(interp,objv[1],&wavelength);//Angstrom
  Tcl_GetDoubleFromObj(interp,objv[2],&abs_length);//mm
  Tcl_GetDoubleFromObj(interp,objv[3],&splength);//not used for flat sample
  Tcl_GetDoubleFromObj(interp,objv[4],&/*convert*/polydisp);//0/1
  //this is used to recognize ULVs and turn on appropriate polydispersity and absorption corrections
  if (polydisp==1) type=4;	//ULVs
  if (polydisp == 2) type = 5; //Input F[h] is already form factor and no correction is necessary
  Tcl_GetDoubleFromObj(interp,objv[5],&thickness);//micron
  Tcl_GetDoubleFromObj(interp,objv[6],&lorentz);//power
  Tcl_GetDoubleFromObj(interp,objv[7],&D);//Angstrom
  Tcl_GetDoubleFromObj(interp,objv[8],&A);//Angstrom squared
  thickness*=0.001;
  hn=(objc-19)/3;
  printf("%d %d********1\n",hn,objc);
  free(F); F=(double*)malloc(sizeof(double)*(hn+1));
  free(Q); Q=(double*)malloc(sizeof(double)*(hn+1));
  free(errf); errf=(double*)malloc(sizeof(double)*(hn+1));
  free(TOT); TOT=(double*)malloc(sizeof(double)*(hn+1));
  Tcl_GetDoubleFromObj(interp,objv[9],&F[0]);

  char* str; str=Tcl_GetString(objv[10]);
  if (strlen(str)>1) contrast=atoi(&str[1]);
  mode=str[0];
  if (mode=='n') {
	  frmelm.graph=ctfgraph[0];
	  frmelm.graphr=rhograph[0];
  } else if (mode=='x') {
	  frmelm.graph=ctfgraph[1];
	  frmelm.graphr=rhograph[1];
  }
  Tcl_GetDoubleFromObj(interp,objv[11],&rw);//Angstrom^-2
  Tcl_GetDoubleFromObj(interp,objv[12],&eCG);//Angstrom or electrons
  Tcl_GetDoubleFromObj(interp,objv[13],&ePh);
  Tcl_GetDoubleFromObj(interp,objv[14],&eCh);
  Tcl_GetDoubleFromObj(interp,objv[15],&ec2);
  Tcl_GetDoubleFromObj(interp,objv[16],&ec1);
  Tcl_GetDoubleFromObj(interp,objv[17],&ec3);
  Tcl_GetDoubleFromObj(interp,objv[18],&epep);

  j=19;
  Q[0] = 0;
  for(i=1;i<=hn;i++){
    Tcl_GetDoubleFromObj(interp,objv[j++],&F[i]);
    Tcl_GetDoubleFromObj(interp,objv[j++],&Q[i]);
    Tcl_GetDoubleFromObj(interp,objv[j++],&errf[i]);
    if(direct_err == 0) errf[i]+=sqrt(fabs(F[i]*stepsize_integral));
  }
	
  //if(D<HNUM && D>0) D=(D)*wavelength/2./sin(Q[int(D)]/2.*PI/180.);
  if(D>0) for(i=1;i<=hn;i++) Q[i]=2.*PI*i/D;
  else D=fabs(D);
  correction();
  //	normalize();
//  free(TOT);TOT=NULL;
  Blt_ResizeVector(frmelm.xv,hn+1);
  Blt_ResizeVector(frmelm.yv,hn+1);

  return(0);
}

double sample::interpolq2f(double q){
  double newF=0, temp;
  int k=1;
  while(k<=hn && Q[k]<q) k++;
  if(k>hn) printf("Error: q is outside the range of the simulated curve\n",q);
  else polint(Q+k-0, TOT+k-0, 1, q, &newF, &temp);
  return newF;
}

int sample::cft(double *qP,double *fP){
  /* prepare for plotting fitting curve */
  double q,step=1.4/2000;
  int i=0;
  for(q=0.0001;q<1.4;q+=step){//1999 pts
    qP[i]=q;
		if(showsign) fP[i++]=(numerical==0?q2f(q):q2fnumer(q));
        else fP[i++]=fabs((numerical==0?q2f(q):q2fnumer(q)));
  }
  return 1999;
}

int sample::sdp(double *zP,double *rP){
  /* use z2rho to construct the entire SDP of the model */
  int i,PT=1000;
  double STEP,z;
  STEP=g2.D/PT/2.0;
//  amplitudes(); make sure it is called before sdp calculation
  for(i=0,z=0;i<=PT;i++,z+=STEP){
    rP[i+PT]=z2sdp(z);
    zP[i+PT]=z;
  }
  for(i=0;i<PT;i++) {
    zP[i]=-zP[2*PT-i];
    rP[i]=rP[2*PT-i];
  }
  return 2001;
}


double sample::z2sdp(double z){
  /* given z return \sdp for one point*/
  double f;
//  amplitudes(); moved to sdp
  f= ((g2.VCG==0?0:eCG/ g2.VCG)-rw) * g2.z2rhoCG(z)
	+((g2.VPh==0?0:ePh/ g2.VPh)-rw) * g2.z2rhoPh(z)
	+((g2.VCh==0?0:eCh/ g2.VCh)-rw) * g2.z2rhoCh(z)
	+((g2.Vc2==0?0:ec2/ g2.Vc2)-rw)* g2.z2rhoc2(z)
	+((g2.Vc1==0?0:ec1/ g2.Vc1)-rw)* g2.z2rhoc1(z)
	+((g2.Vc3==0?0:ec3/ g2.Vc3)-rw)* g2.z2rhoc3(z)
	+((g2.Vpep==0?0:epep/ g2.Vpep)-rw) * g2.z2rhopep(z)
	+rw;
  return f;
}


double sample::q2f(double q){
  // given q return F(q)
  double XCGq,XPhq,XChq,XCq,Xc1q,Xc3q;
  double SIGCGq,SIGPhq,SIGChq,SIGCq,SIGc1q,SIGc3q;
//  amplitudes(); make sure that amplitudes are called before q2f calculation
  XCGq=g2.parCG[0]*q;
  XPhq=g2.parPh[0]*q;
  XChq=g2.parCh[0]*q;
  XCq=g2.parC[0]*q;
  Xc1q=g2.parc1[0]*q;
  Xc3q=g2.parc3[0]*q;
  SIGCGq=g2.parCG[2]*q;
  SIGPhq=g2.parPh[2]*q;
  SIGChq=g2.parCh[2]*q;
  SIGCq=g2.parC[2]*q;
  SIGc1q=g2.parc1[2]*q;
  SIGc3q=g2.parc3[2]*q;

  return( 2*( ((g2.VCG==0?0:eCG/ g2.VCG)-rw) *g2.parCG[1] * g2.parCG[2]*cos(XCGq)*exp(-0.5*SIGCGq*SIGCGq)
		  +   ((g2.VPh==0?0:ePh/ g2.VPh)-rw) *g2.parPh[1] * g2.parPh[2]*cos(XPhq)*exp(-0.5*SIGPhq*SIGPhq)
		  +   ((g2.VCh==0?0:eCh/ g2.VCh)-rw) *g2.parCh[1] * g2.parCh[2]*cos(XChq)*exp(-0.5*SIGChq*SIGChq)
		  +   ((g2.Vc3==0?0:ec3/ g2.Vc3)-rw) *g2.parc3[1] * g2.parc3[2]*cos(Xc3q)*exp(-0.5*SIGc3q*SIGc3q) 
		  +   ((g2.Vc1==0?0:ec1/ g2.Vc1)-rw) *g2.parc1[1] * g2.parc1[2]*cos(Xc1q)*exp(-0.5*SIGc1q*SIGc1q) 
		  + ( ((g2.Vc2==0?0:ec2/ g2.Vc2)-rw) * g2.parC[1] + ((g2.Vpep==0?0:epep/g2.Vpep )-rw)*(1-g2.parC[1]) )
			 * ((sin(XCq)*exp(-0.5*SIGCq*SIGCq)/q ) - g2.parc3[1]*g2.parc3[2]*cos(Xc3q)*exp(-0.5*SIGc3q*SIGc3q)
													- g2.parc1[1]*g2.parc1[2]*cos(Xc1q)*exp(-0.5*SIGc1q*SIGc1q)) ) );
}


double sample::q2fnumer(double q){
  // given q return F(q) calculated numerically
    int j,PT=1000;
    double Re=0,STEP=D/PT/2.0,z;
    for (j=0,z=0;j<=PT;j++,z+=STEP) {
      Re+=(z2sdp(z)-rw)*cos(q*z);
    }
    return 2*Re*STEP;  
}


int sample::assignsigns(int objc, Tcl_Obj *const objv[]){
  /* this assigns the signs to the forfactor values */
  int sign;
  for (int i=0; i<objc; i++){
      Tcl_GetIntFromObj(interp,objv[i],&sign);
	  F[i]=sign*fabs(F[i]);
  }
  for(int h=1;h<=hn;h++){ TOT[h]=F[h]; }
  return(0);
}


int sample::RHOdft(double *zP,double *rP){
  /* calculate profile using discrete Fourier transform */
  double z,dz,hg=0,*qs,*fs;
  int i,h;
  const int PT=1000;
  double Dtemp=D>100?g2.D:D;
  qs=(frmelm.xv)->valueArr;
  fs=(frmelm.yv)->valueArr;
  for(i=0,z=0,dz=Dtemp/PT/2;i<=PT;i++,z+=dz){
    for(h=1,rP[i+PT]=0;h<=hn;h++) rP[i+PT]+=F[h]*cos(Q[h]*z);
    rP[i+PT]=(rP[i+PT]*2*SCL+F[0])/D;
    zP[i+PT]=z;
//    if(rP[i+PT]>0) hg+=rP[i+PT];
  }
//  printf("%lf\n",hg*A*dz);
  for(i=0;i<PT;i++) {
    zP[i]=-zP[2*PT-i];
    rP[i]=rP[2*PT-i];
  }
  return 2001;
}


int sample::dft(double *qP, double *fP, double *zP, double *rP){
  /* prepare for plotting fitting curve */
  double q,step=1.4/2000;
  int i=0;
  printf("SCL = %lf\n",SCL);
  for(q=0.0001;q<1.4;q+=step){//1999 pts
    qP[i]=q;
		if(showsign) fP[i++]=(q2df(q,zP,rP));
        else fP[i++]=fabs(q2df(q,zP,rP));
  }
  return 1999;
}


double sample::q2df(double q, double *zP, double *rP){
  // given q return F(q) calculated numericaly from dft profile
    int j,PT=1000;
    double Re=0,STEP=D/PT/2.0;
    for (j=PT;j<=PT+PT;j++) {
      Re+=(rP[j]-rP[0])*cos(q*zP[j]);
    }
    return 2*Re*STEP;  
}


int sample::frm(double *qP,double *fP){
  // prepare for plotting data form factors 
  for(int i=0;i<=hn;i++){
//    NK's modification for showing the negative exp. values
//    if(showsign) fP[i]=(F[i]);
//    else fP[i]=fabs(F[i]);
    fP[i]=(F[i]);
    qP[i]=Q[i];
  }
  return hn+1;
}


void sample::drawfrm(){
  /* draw data form factors of the sample */
  double *qs=(frmelm.xv)->valueArr, *fs=(frmelm.yv)->valueArr;
  int points;
//  amplitudes();
  if(frmelm.status==0){
    points=frm(qs,fs);
    if(hn>60){
      Blt_ResetVector(frmelm.xv,qs, hn+1, hn+1, TCL_VOLATILE);
      Blt_ResetVector(frmelm.yv,fs, hn+1, hn+1, TCL_VOLATILE);
    }else{
      Blt_ResetVector(frmelm.xv,qs, hn+1, 32, TCL_STATIC);
      Blt_ResetVector(frmelm.yv,fs, hn+1, 32, TCL_STATIC);
    }
    sprintf(cmd,"%s element create %d -pixels 2 -linewidth 0 -xdata xfrm%d -ydata yfrm%d -color %s",frmelm.graph,frmelm.id,frmelm.id,frmelm.id,colors[frmelm.id]);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);

	qs=(frmelm.xmv)->valueArr;
	fs=(frmelm.ymv)->valueArr;
//	points=cft(qs,fs);
	Blt_ResetVector(frmelm.xmv,qs, 1, 1999, TCL_VOLATILE);
	Blt_ResetVector(frmelm.ymv,fs, 1, 1999, TCL_VOLATILE);
    sprintf(cmd,"%s element create m%d -pixels 0 -linewidth 2 -xdata xmdc%d -ydata ymdc%d -color %s",frmelm.graph,frmelm.id,frmelm.id,frmelm.id,colors[frmelm.id]);
//    sprintf(cmd,"%s element create m%d -pixels 0 -linewidth 2 -xdata xmdc%d -ydata ymdc%d -color black",frmelm.graph,frmelm.id,frmelm.id,frmelm.id);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);

	qs=(frmelm.xmsp)->valueArr;
	fs=(frmelm.ymsp)->valueArr;
	points=sdp(qs,fs);
	Blt_ResetVector(frmelm.xmsp,qs, points, points, TCL_VOLATILE);
	Blt_ResetVector(frmelm.ymsp,fs, points, points, TCL_VOLATILE);
    sprintf(cmd,"%s element create s%d -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymsp%d -color %s",frmelm.graphr,frmelm.id,frmelm.id,frmelm.id,colors[frmelm.id]);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);

    frmelm.status=1;
  // Draw EDP
   //sprintf(cmd,"calculateEDP");
   //Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
   sprintf(cmd,"%s element create CGEDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymCG_EDP_v -color black","$rhowX",frmelm.id);
Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
     sprintf(cmd,"%s element create PhEDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymPh_EDP_v -color red","$rhowX",frmelm.id);
Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
     sprintf(cmd,"%s element create ChEDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymCh_EDP_v -color green","$rhowX",frmelm.id);
Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
	   sprintf(cmd,"%s element create c2EDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymc2_EDP_v -color blue","$rhowX",frmelm.id);
Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
      sprintf(cmd,"%s element create c1EDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymc1_EDP_v -color cyan","$rhowX",frmelm.id); 
Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);  
        sprintf(cmd,"%s element create c3EDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymc3_EDP_v -color magenta","$rhowX",frmelm.id);
   Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
sprintf(cmd,"%s element create wEDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ymw_EDP_v -color yellow","$rhowX",frmelm.id);
   Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
sprintf(cmd,"%s element create pepEDP -pixels 0 -linewidth 2 -xdata xmsp%d -ydata ympep_EDP_v -color #888800","$rhowX",frmelm.id);
   Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);

    sprintf(cmd,"raiseq %s %d",frmelm.graph,frmelm.id);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);    
  }else{
    sprintf(cmd,"%s element delete %d",frmelm.graph,frmelm.id);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    sprintf(cmd,"%s element delete m%d",frmelm.graph,frmelm.id);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
	sprintf(cmd,"%s element delete s%d",frmelm.graphr,frmelm.id);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
	sprintf(cmd,"%s element delete CGEDP",frmelm.graphr);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    sprintf(cmd,"%s element delete PhEDP",frmelm.graphr);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    sprintf(cmd,"%s element delete ChEDP",frmelm.graphr);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    sprintf(cmd,"%s element delete c2EDP",frmelm.graphr);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    sprintf(cmd,"%s element delete c1EDP",frmelm.graphr);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    sprintf(cmd,"%s element delete c3EDP",frmelm.graphr);
    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    sprintf(cmd,"%s element delete wEDP",frmelm.graphr);

    Tcl_EvalEx(interp,cmd,-1,TCL_EVAL_GLOBAL);
    frmelm.status=0;
  }
}


/*********************
some helper functions for different corrections 
************************/

double sample::absorp1(double s)
{
  double ct,path;
  ct=sqrt(1-tmp*tmp);
  path = sqrt(pow(splength+s,2.) - pow(splength*ct,2.));
  return exp( -2*path/abs_length);
}
  
  
sample *dirty;
double absorp1_wrap(double s){
  return dirty->absorp1(s);
}


void sample::curvecor() {
  dirty=this;
  for (int i=1; i<=hn; i++){
    tmp=Q[i]*wavelength/(4*PI);
    TOT[i]=(thickness/exp(2*splength*tmp/abs_length))/qromb(absorp1_wrap,0,thickness)*pow(Q[i],lorentz);   // wiener absorp2
    printf("cur: %f\n", TOT[i]);  
  }
}


void sample::flatcor(){
  /* flat sample corrections
   includes absorption, footprint, Lorentz  corrections */
  FILE *fp;
  double *FP,*FPC,*ABS,xless[2000],inteless[2000];
//  FP=(double*)malloc(sizeof(double)*hn);
//  FPC=(double*)malloc(sizeof(double)*hn);
//  ABS=(double*)malloc(sizeof(double)*hn);
  double allbeam=0,partbeam=0.,temp=1;
  int i=0,j=0,n=0,n1=0,n2=0,n3=0;
  printf("sample # %d>>>>",number);
  fp=fopen(beamfile,"r");
  if(fp!=NULL){
    while(fscanf(fp,"%lf %lf",&xless[j],&inteless[j])!=EOF) {j++;n++;}
  }else {
    printf("error in openning beamfile: %s\n",beamfile);
    for(i=1;i<=hn;i++){
      double temp, temp2, abs;
      temp=Q[i]*wavelength/4/PI;
	  //printf("abs: %f\n",abs_length);
	  if (abs_length<=0) TOT[i]=pow(Q[i],lorentz);
	  else {
		  // Yufeng's formula
		  /*
		  if((temp2=2*thickness/abs_length/temp)>49) temp2=1; else temp2=1-exp(-temp2);
		  abs=(1-exp(-2.e-7*D*50/abs_length/temp))/temp2;
		  TOT[i]=pow(Q[i],lorentz)*abs;
		  */
		  
		  // NK's formula
//		  /*
		  double cons = (normal_mode==0 ? 1e-7*D*50/thickness : 1); // for consistency with the Yufeng's weird (non)normalization
		  temp2=2*thickness/abs_length/temp;
		  abs=cons*temp2/(1-exp(-temp2));
		  TOT[i]=pow(Q[i],lorentz)*abs;
//		  */
	  }
    }
//    free(FP);
//    free(FPC);
//    free(ABS);
    return;
  }
  fclose(fp);
//  free(FP);
//  free(FPC);
//  free(ABS);
}


void sample::ulvcor(){
  /* ulv sample corrections
   includes absorption, Lorentz  corrections */
  printf("sample # %d>>>>",number);
  for(int i=1;i<=hn;i++){
	  double temp, temp2, abs;
      temp=Q[i]*wavelength/4/PI;
	  temp=1-2*temp*temp;
	  if (abs_length<=0) TOT[i]=pow(Q[i],lorentz);
	  else {
		  // NK's formula
//		  /*
		  temp2=(1-temp)/temp/abs_length*thickness;
		  abs=temp2/(1-exp(-temp2));///exp(-thickness/abs_length);	dropped due to the normalization purposes
		  TOT[i]=pow(Q[i],lorentz)*abs;
//		  */
	  }
  }
}


void sample::correction(){
  double temp=1.0;
  int h;
  switch(type){
  case 1:
    flatcor();
	break;
  case 2:
    curvecor();
	break;
  case 4:
    ulvcor();
	break;
  case 5:
    //In this case, the input smp file contains form factor, not scaling factor.
    for (h = 1; h <= hn; h++) {
      TOT[h] = F[h];
    }
  }
  for (h=1;h<=hn;h++){
    F[h]=(F[h]<0?-1:1)*sqrt(fabs(F[h]*TOT[h]));
    if(direct_err == 0) errf[h]=(sqrt(errf[h]*TOT[h]+F[h]*F[h])-fabs(F[h]));
  }
  cout << "TOT[1] = " << TOT[1] << endl;
}

/*  This ends the corrections  */
