#include "SysDefs.h"

void runge(n,h,x,t, sysdata)
	int n;		/* number of the equations */
	double h;	/* width */ 
	double x[];	/* input/output coordinate with initial value */
	double t;	/* time */
	SysData *sysdata;
{
	double k1[NN],k2[NN],k3[NN],k4[NN];
	double xtemp[NN];
	int i;

	function(x,k1,t, sysdata);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k1[i]/2.0*h;
	t+=h/2.0;

	function(xtemp,k2,t, sysdata);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k2[i]/2.0*h;

	function(xtemp,k3,t, sysdata);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k3[i]*h;
	t+=h/2.0;

	function(xtemp,k4,t, sysdata);
	for (i=0; i<n; i++) x[i]+=(k1[i]+2.0*(k2[i]+k3[i])+k4[i])/6.0*h;

}
