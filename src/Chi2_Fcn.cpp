#include "Chi2_Fcn.h"
#include "g2model.h"
#include "

double Chi2_Fcn::operator()()
{  
  int i,j;

  for(i=0;i<rdim;i++) x[xindex[i]]=pars[i];
  x[0]=fabs(x[0]);x[2]=fabs(x[2]);
  x[3]=fabs(x[3]);x[5]=fabs(x[5]);
  x[6]=fabs(x[6]);x[8]=fabs(x[8]);
  x[9]=fabs(x[9]);x[11]=fabs(x[11]);
  x[12]=fabs(x[12]);x[14]=fabs(x[14]);
  x[15]=fabs(x[15]);x[16]=fabs(x[16]);

  if (x[4]*x[5]==0) {x[3]=x[0];}

  g2.penalty = penalty();
  g2.basekN = 0, g2.basekX = 0;
  g2.rsdtN = 0, g2.rsdtX = 0;
  for (j = 0; j < SNUM; j++) {
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
