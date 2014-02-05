#ifndef GUARD_G2MODEL_H
#define GUARD_G2MODEL_H


class Gauss {
public:
  void set_params(double a, double b, double c) { 
    center = a; ampl = b; sigma = c; 
  }
  void set_targets(double a, double b, double c) {
    target_c = a; target_a = b; target_s = c;
  }                
  void set_tols(double a, double b, double c) {
    tol_c = a; tol_a = b; tol_s = c;
  }
  double penalty();
  double center; // center
  double ampl; // amplitude
  double sigma; // sigma
  double target_c;
  double target_a;
  double target_s;
  double tol_c;
  double tol_a;
  double tol_s; 
};


class G2model {
// delaration of the data structure representing the model of a bilayer profile
// in terms of (volume) probability distribution functions
protected:
public:
  char name[16];
  double parCG[3]; //XH1,C1,SIG1			- CG gaussian			//0-2
  double parPh[3]; //XH2,C2,SIG2			- PCN gaussian			//3-5
  double parCh[3]; //XH3,C3,SIG3			- ChCH3 gaussian		//6-8
  double parC[3]; //XC,CC,SIGC			  - C error function		//9-11
  double parc1[3]; //XDB,CDB,SIGDB		- double bond gaussian	//12-14
  double parc3[3]; //XM,CM,SIGM			- methyl gaussian		//15-17
  double parr,parr12; // r = Vc3/Vc2, r12 = Vc1/Vc2				//18,19
  double RCG,RPh; // RCG = VCG/VHL, RPh = VPh/VHL				//20,21
  double Rm,sigR; //polydispersity								//22,23

  double A;														//24

  double /*e1,e2,*/nc2,nc1,nc3/*,ec2,em,rw,ep*/;				//25-27
  double VL,VHL,Vpep;											//28-30
  double VH,VCG,VPh,VCh,Vc2,Vc1,Vc3;							//31-37
  double D,DHH,DB,DC;											//38-41
  double DH1,tDH1,DH1R;											//42-44
  double dXH,tdXH,sdXH;											//45-47
  double tr,sr,tr12,sr12;										//48-51
  double tRCG,sRCG,tRPh,sRPh;									//52-55
  double tSIGC,sSIGC;											//56,57
  double tDC,sDC;												//58,59
  double dXH2,tdXH2,sdXH2;										//60-62

  double basek,rsd;												//63,64
  double wiggl,twiggl;											//65,66
  double dXH3,tdXH3,sdXH3;										//67-69

  double rsdj,basekN,basekX;
  double rsdtN,rsdtX;
  double rhoc2;
  int hn;
  int rho(double *z,/*double *,*/double *pCG,double *pPh,double *pCh,double *pc2,double *pc1,double *pc3,double *pW,double *ppep);
  void init();
  void linkvar();
  double max();
//  double wiggeling(double *);
  double wiggeling();
//  double gibbs();
  double z2rhoCG(double);
  double z2rhoPh(double);
  double z2rhoCh(double);
  double z2rhoc1(double);
  double z2rhoc2(double);
  double z2rhoc3(double);
  double z2rhoC(double);
  double z2rhoW(double);
  double z2rhopep(double);
  double penalty;
  Gauss carbGlyc; // CG = Carbonyl-Glycerol group
  Gauss phosphate; // Ph
  Gauss choline; // Ch
  Gauss methine; // c1 = double bond, CH
  Gauss methyl; // c3 = terminal CH3
  double get_penalty();
};

#endif
