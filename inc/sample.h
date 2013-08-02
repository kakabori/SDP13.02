#ifndef GUARD_sample_h
#define GUARD_sample_h

#include <blt.h>
#include <stdlib.h>
#include <iostream>
#include "bilayer.h"

class sample{
  /* data structure represents a sample data set */
public:
  sample();
  ~sample();
  double tmp, D, P, A, SCL;
  double wavelength, lorentz, abs_length, splength, polydisp, thickness;
  char mode;
  int contrast;
  double rw, eCG, ePh, eCh, ec1, ec2, ec3, epep;
  int hn, type;
  double basek;
  //*F is an array to store Lorentz and absorption corrected form factor
  double *F,*Q,*errf, *TOT;
  //*Phi is an array to store the original scaling factor
  //*sigmaPhi is an array to store the original uncertainties on the scaling factor
  //*absLorCorr is an array to store the absorption and Lorentz correction combined factor
  double *Phi, *sigmaPhi, *absLorCorr;
  char color[16];
  char name[32];
  char *beamfile;
  int number;
  Bltelm frmelm;
  int assignvalues(int objc, Tcl_Obj *const objv[]);
  double interpolq2f(double q);
  void initialize(char *);
  int frm(double *qP,double *fP);
  int cft(double *qP,double *fP);
  int sdp(double *qP,double *fP);
  int assignsigns(int objc, Tcl_Obj *const objv[]);
  int RHOdft(double *qP,double *fP);
  int dft(double *qP,double *fP, double *zP, double *rP);
  double q2df(double q, double *zP, double *rP);
  void drawfrm();
  void correction();
  void flatcor();
  void ulvcor();
  void curvecor();
  void normalize();
  void scale2model();
  void scale2sim();
  void redrawfrm();
  double z2sdp(double);
  double q2f(double);
  double q2fnumer(double);
  double absorp2(double s);
  double absorp1(double s);
  double absconst(double s);
private:
};

#endif
