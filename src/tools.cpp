
/*********************************************

utility functions from Numerical Recipe

includes polynomial interpolation
and numerical integration (qromb)

*********************************************/
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <malloc.h>

#define DF_PREC double
#define Pi 3.1415926535897932384

void nrerror(char *error_text);
DF_PREC *vector(int nl,int nh);
void free_vector(DF_PREC *v,int nl,int nh);

static DF_PREC qromb_ss,qromb_dss;
static DF_PREC qromb_s[32],qromb_h[32]={1,1,1};
static int qromb_j;

void nrerror(const char *error_text){
  fprintf(stderr,"%s\n",error_text);
}


DF_PREC *vector(int nl,int nh){
  DF_PREC *v;
  v=(DF_PREC *)malloc((unsigned) (nh-nl+1)*sizeof(DF_PREC));
  if (!v) nrerror("allocation failure in vector ()");
  return v-nl;
}


void free_vector(DF_PREC *v,int nl,int nh){
  free((char*) (v+nl));
}


DF_PREC polint4s(DF_PREC xmx0,double *ya){
  return (ya[0]*(6-xmx0)+ya[3]*xmx0)*(xmx0-2)*(xmx0-4)/48+
    (ya[2]*(2-xmx0)+ya[1]*(xmx0-4))*(xmx0-6)*xmx0/16;
}


DF_PREC polint4t(DF_PREC xmx0,double *ya){
  return (ya[0]*(3-xmx0)+ya[3]*xmx0)*(xmx0-1)*(xmx0-2)/6+
    (ya[2]*(1-xmx0)+ya[1]*(xmx0-2))*(xmx0-3)*xmx0/2;
}


DF_PREC polint3(int xmx1,int xmx2,int xmx3,double y1,double y2,double y3){
  DF_PREC x2mx3,x3mx1;
  x2mx3=xmx3-xmx2;
  x3mx1=xmx1-xmx3;
  return (xmx2*x2mx3*(xmx3*y1-xmx1*y3)+xmx1*x3mx1*(xmx3*y2-xmx2*y3))/(x2mx3*x3mx1*(x2mx3+x3mx1));
}


DF_PREC polint3_log(DF_PREC xmxa,float *ya){
  return xmxa*50*((100*xmxa-2)*(ya[0]+ya[2]-2*ya[1])+ya[2]-ya[0])+ya[0];
}


DF_PREC polint30_log(DF_PREC x,float *ya){
  return (0.9772372212*ya[0]*(x-1)-42.93136751*ya[1]*x)*(x-1.023292992)+41.95413029*ya[2]*(x-1)*x;
}


void polint(DF_PREC *xa,DF_PREC *ya,int n,DF_PREC x,DF_PREC *y,DF_PREC *dy)
{
  int i,m,ns=1;
  DF_PREC den,dif,dift,ho,hp,w;
  DF_PREC *c,*d;
  dif=fabs(x-xa[1]);
  c=vector(1,n);
  d=vector(1,n);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) nrerror("Error in routine POLINT");
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_vector(d,1,n);
  free_vector(c,1,n);
}


#define FUNC(x) ((*func)(x))

DF_PREC trapzd(DF_PREC (*func)(DF_PREC),DF_PREC a,DF_PREC b,int n)
{
  register DF_PREC x,sum,del;
  DF_PREC tnm;
  static DF_PREC s;
  static int it;
  int j;
  
  if (n == 1) {
    it=1;
    return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
  } else {
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
    it *= 2;
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}


//#define EPS 1.0e-6
#define JMAX 22
#define JMAXP JMAX+1
#define K 11
#define KT 5
DF_PREC qromb(DF_PREC (*func)(DF_PREC),DF_PREC a,DF_PREC b)
{
  double EPS=1.0e-6;
  for (qromb_j=1;qromb_j<=JMAX;qromb_j++) {
    qromb_s[qromb_j]=trapzd(func,a,b,qromb_j);
    if (qromb_j >= K) {
      polint(&qromb_h[qromb_j-K+KT],&qromb_s[qromb_j-K+KT],K-KT,0.0,&qromb_ss,&qromb_dss);
      //			cout<<qromb_ss<<' '<<qromb_dss<<endl;
      if (fabs(qromb_dss) < EPS*fabs(qromb_ss)) {
	//				cout<<"loops   "<<qromb_j<<endl;
	return qromb_ss;
      }
    }
    qromb_h[qromb_j+1]=0.25*qromb_h[qromb_j];
  }
  nrerror("Too many steps in routine QROMB");
  return qromb_ss;
}


//************************************//
// error function calculation from NR //
//************************************//
#define ITMAX 100 //Maximum allowed number of iterations.
//#define EPS 3.0e-7 //Relative accuracy.
#define FPMIN 1.0e-30 //Number near the smallest representable floating-point number.
double gammln(double xx)
//Returns the value ln[?(xx)] for xx > 0.
{
//Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
//accuracy is good enough.
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


void gser(double *gamser, double a, double x, double *gln)
//Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
//Also returns ln ?(a) as gln.
{
//	double gammln(double xx);
	int n;
	double sum,del,ap;
	double EPS=3.0e-7;
	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) 	printf("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf("a too large, ITMAX too small in routine gser");
		return;
	}
}


void gcf(double *gammcf, double a, double x, double *gln)
//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
//as gammcf. Also returns ln?(a) as gln.
{
//	double gammln(double xx);
	int i;
	double an,b,c,d,del,h;
	double EPS=3.0e-7;
	*gln=gammln(a);
	b=x+1.0-a; //Set up for evaluating continued fraction by modified Lentz’s method (§5.2)	with b0 = 0.
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) { //Iterate to convergence.
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) printf("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h; //Put factors in front.
}


double gammp(double a, double x)
//Returns the incomplete gamma function P(a, x).
{
//	void gcf(double *gammcf, double a, double x, double *gln);
//	void gser(double *gamser, double a, double x, double *gln);
	double gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) printf("Invalid arguments in routine gammp");
	if (x < (a+1.0)) { //Use the series representation.
		gser(&gamser,a,x,&gln);
		return gamser;
	} else { //Use the continued fraction representation.
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf; //and take its complement.
	}
}


double erff(double x)
//Returns the error function erf(x).
{
//	double gammp(double a, double x);
	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}


#undef ITMAX
//#undef EPS
#undef FPMIN
//************************************//


double NKerff(double x)
//Returns the error function value from table
{
  int sign=1;
  double NKerr[200]={0,0.01701,0.03401,0.051,0.06796,0.08489,0.10179,0.11863,0.13543,0.15216,
  0.16883,0.18542,0.20192,0.21834,0.23466,0.25088,0.26698,0.28297,0.29884,0.31458,
  0.33018,0.34564,0.36096,0.37612,0.39112,0.40596,0.42064,0.43514,0.44946,0.46361,
  0.47756,0.49133,0.50491,0.51829,0.53147,0.54445,0.55722,0.56979,0.58215,0.59429,
  0.60623,0.61794,0.62944,0.64073,0.65179,0.66264,0.67326,0.68367,0.69386,0.70382,
  0.71357,0.7231,0.73241,0.7415,0.75038,0.75904,0.76749,0.77572,0.78375,0.79156,
  0.79917,0.80657,0.81377,0.82078,0.82758,0.83419,0.8406,0.84683,0.85287,0.85873,
  0.8644,0.8699,0.87522,0.88037,0.88536,0.89018,0.89483,0.89933,0.90368,0.90787,
  0.91191,0.91582,0.91957,0.9232,0.92668,0.93004,0.93327,0.93638,0.93936,0.94223,
  0.94499,0.94763,0.95017,0.9526,0.95494,0.95717,0.95931,0.96136,0.96332,0.9652,
  0.96699,0.9687,0.97034,0.9719,0.97339,0.97482,0.97617,0.97746,0.9787,0.97987,
  0.98098,0.98204,0.98305,0.98401,0.98492,0.98578,0.98661,0.98738,0.98812,0.98882,
  0.98948,0.99011,0.99071,0.99127,0.9918,0.9923,0.99278,0.99322,0.99365,0.99405,
  0.99442,0.99478,0.99511,0.99543,0.99572,0.996,0.99626,0.99651,0.99674,0.99696,
  0.99716,0.99735,0.99753,0.9977,0.99786,0.99801,0.99815,0.99828,0.9984,0.99851,
  0.99862,0.99871,0.99881,0.99889,0.99897,0.99905,0.99912,0.99918,0.99924,0.9993,
  0.99935,0.9994,0.99945,0.99949,0.99953,0.99956,0.9996,0.99963,0.99966,0.99969,
  0.99971,0.99973,0.99975,0.99977,0.99979,0.99981,0.99982,0.99984,0.99985,0.99986,
  0.99988,0.99989,0.9999,0.9999,0.99991,0.99992,0.99993,0.99993,0.99994,0.99994,
  0.99995,0.99995,0.99996,0.99996,0.99996,0.99997,0.99997,0.99997,0.99998,0.99998};
  
  if (x<0.0) {x*=-1.0; sign=-1;}
  double x0=x*66;
  if (x>3.0) return 1*sign;
  else return sign*(NKerr[int(x0)]+(NKerr[int(x0)+1]-NKerr[int(x0)])*(x0-int(x0)));
}
