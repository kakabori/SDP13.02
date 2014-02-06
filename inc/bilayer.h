#ifndef GUARD_BILAYER_H
#define GUARD_BILAYER_H

double int_to_ext(unsigned int, double);
double ext_to_int(unsigned int, double);
int setNePep(ClientData, Tcl_Interp *, int, Tcl_Obj *const objv[]);
double amotry(double **p, double *y, double *psum, int ndim, double (*funk)(double *), int ihi, double fac);
int amoeba(double **p, double *y, int ndim, double ftol, double (*funk)(double *));
int MR_calculateEDP(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]);
void after_fitting();
double ran1(long *idum);
double active_ctf(double);
double active_interpol(double);
double qromb(double (*Func)(double), double a, double b);
double erff(double x);
double NKerff(double x);
double polydispcorr(double, double, double);
void polint(const double *xa, const double *ya,int n,double x,double *y,double *dy);

struct Bltelm {
  int id;
  char *graph;
  char *graphr;
  char status;
  Blt_Vector *xv, *yv, *xmv, *ymv;//experimental FormFactors and Model fit
  Blt_Vector *xmsp, *ymsp;//models of scattering density profiles
};

#endif
