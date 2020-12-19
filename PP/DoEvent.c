#include "SysDefs.h"

void DoEvent(double x[], double t, SysData *sys)
{
	int i, test, diff;

	switch (xddpret.key){
		case ' ':	
			xhplot(CLS, x[0],x[1],BLACK);
			xhplot(FRAME, x[0],x[1],WHITE);
		break;	
		case 's':	/* letter s */
			printf("Status\n");
			for (i = 0; i < NP; i++){
				printf("P[%d] = %lf, ",i, sys->param[i]);
			}
			printf("\n");
			printf("location( ");
			for (i = 0; i < NN; i++){ printf("%lf ", x[i]); }
			printf(")\n");
			printf("Status\n");
			for (i = 0; i < NP; i++){ printf("%.10f ", sys->param[i]); }
			for (i = 0; i < NN; i++){ printf("%.15f ", x[i]); }
			printf("\n");
		break;
		case 'f':	/* letter f */
			sys->pflag = 1 - sys->pflag;
			printf("Poincare sw\n");
		break;
		case 'q':	/* letter q */
			printf("quit\n");
			exit(0);
		break;
		case 'i':
			x[0] = sys->x0[0][0];
			x[1] = sys->x0[0][1];
			x[2] = sys->x0[0][2];
		break;
	}
	test = xddpret.key - 'A'; 
	diff = 'a' - 'A';

	if (test >= 0 && test < NP){
		sys->param[test] += sys->dparam[test];
		printf("param[%d] = %lf \n", test, sys->param[test]);
	}
	if (test >= diff && test < NP + diff){
		sys->param[test-diff] -= sys->dparam[test-diff];
		printf("param[%d] = %lf \n", test-diff, sys->param[test-diff]);
	}
	return ;
}

