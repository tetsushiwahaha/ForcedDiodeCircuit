#ifdef __FreeBSD__
#include <floatingpoint.h>
#endif

#include "SysDefs.h"

int main(argc, argv)
	int argc; 
	char **argv;
{
	int i, j, k;
	int ii = 0;
	double x[NN], xprev[NN];
	double t = 0.0, time, tend, h;
	SysData *dpp; 
	int loop;

	if ((dpp = StoreSysData(argc, argv))== NULL){
		exit(-1);
	}

	xhplot(DIRECT_RV, 0.0, 0.0, WHITE);
	
	ii = 0;
	tend = dpp->tend[ii];
	for (i = 0; i < NN; i++){ xprev[i] = x[i] = dpp->x0[ii][i];}
	h = dpp->h;

	xddp.line_wid = dpp->line_wid[0];
	xhplot(PSET, x[0],x[1], WHITE);
	xhplot(LINEATT);

	h = 2.0 * M_PI / dpp->m;
	h *= dpp->tsign[ii];
	while(1){
		for (loop = 0; loop < dpp->m; loop++){
			runge(NN, h, x, t, dpp);
			t += h;
			if (dpp->pflag) {
				xhplot(LINE,x[0],x[1],WHITE);
			}
			if (dpp->write && dpp->pflag){
				fprintf(dpp->write_fp,"%lf %lf %lf\n",x[0],x[1],x[2]);
				fflush(dpp->write_fp); 
			} 
		}

		xhplot(KEYIN);
		if (xddpret.key != -1){ DoEvent(x, t, dpp); }

		if (dpp->write && dpp->pflag == 0){
			fprintf(dpp->write_fp,"%lf %lf\n",x[0],x[1]);
			fflush(dpp->write_fp); 
		} 
		xhplot(ARC,x[0],x[1],WHITE);

		xhplot(POINT, x[0],x[1]);
		if (!xddpret.key){
			x[0] = xddpret.x; x[1] = xddpret.y; 
		}	
		xhplot(PSET, x[0],x[1], BLACK);
	}
}


void function(double x[], double f[], double t, SysData *dpp)
{

	/* positive time simulation. */

	f[0] = x[1] + exp(dpp->param[GAMMA]* x[0]);
	f[1] =  -x[0] - dpp->param[K] * x[1] + dpp->param[B0] 
		+ dpp->param[B] * cos(t);
}

