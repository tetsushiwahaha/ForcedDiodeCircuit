#include <stdio.h>
#include <math.h>
#define NN 6
#define NNN 16

struct forcedDstruct { 
	double gamma;
	double dgamma;

	double k;
	double dk;

	/* direct current term (B0) */
	double B0;
	double dB0;

	/* external alternate (B) */
	double B;
	double dB;

	double eps;			/* eps for newton's method */
	double emax;		/* error maximum boundary */

	int m;				/* integral strip */
	int l;				/* order of period */
	int extra;			/* extrapolation */
	int ite_max;		/* iteration maximum */

	double x0[NN];		/* initial value */
	double x[NN];		/* state variables */
	double xn[NNN];		/* periodic solution */
	double yn[NNN];		/* periodic solution */
	double z[4];		/* characteristic variable */

} forcedD;



main(argc,argv)
	int argc;
	char **argv;
{
	
	int i, j, k = 0;
	int iteration;
	int extra = 0;
	FILE *fopen(), *fpin;
	char c;
	double x_old[2];
	char *tmpname;
	char cmd[BUFSIZ];
	char out[BUFSIZ];

	if (argc == 1){
		fprintf(stderr, "usage: %s filename\n", argv[0]);
		exit(0);
	}

	sprintf(cmd, "/usr/bin/cpp -P %s", argv[1]);

	if ((fpin = popen(cmd,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(-1);
	}

	unlink(tmpname);

	fscanf(fpin,"%lf",&forcedD.gamma);
	fscanf(fpin,"%lf",&forcedD.k);
	fscanf(fpin,"%lf",&forcedD.B0);
	fscanf(fpin,"%lf",&forcedD.B);

	for (i = 0; i < 2; i++){ fscanf(fpin,"%lf",&forcedD.x0[i]); }

	fscanf(fpin,"%lf",&forcedD.dgamma);
	fscanf(fpin,"%lf",&forcedD.dk);
	fscanf(fpin,"%lf",&forcedD.dB0);
	fscanf(fpin,"%lf",&forcedD.dB);

	fscanf(fpin,"%lf",&forcedD.eps);
	fscanf(fpin,"%lf",&forcedD.emax);
	fscanf(fpin,"%d",&forcedD.l);
	fscanf(fpin,"%d",&forcedD.m);
	fscanf(fpin,"%d",&forcedD.ite_max);
	fscanf(fpin,"%d",&forcedD.extra);
	for (i = 0; i < 2; i++){ forcedD.x[i] = forcedD.x0[i]; }

	printf("number of period = %d\n", forcedD.l);

	printf("i gamma    k      b0      b        x(1)      x(2)      mu1     mu2     sqrt\n"); 
	
	while (1) {
		iteration =	fixed(forcedD);
		if (forcedD.z[2] == 0.0){ c = 'C'; } else { c = 'R'; }
		printf("%d %.5f %.5f %.5f %.5f %.7f %.7f %.5f %.5f %.5f %1c\n",
			iteration, forcedD.gamma, forcedD.k, forcedD.B0, forcedD.B, 
				forcedD.x[0], forcedD.x[1], forcedD.z[0], forcedD.z[1],
				sqrt(forcedD.z[0]*forcedD.z[0] + forcedD.z[1]*forcedD.z[1]), 
				c);
		if ( forcedD.l != 1){
			for (j = 0; j < forcedD.l; j++){
				printf("     %8.5f %8.5f\n", forcedD.xn[j], forcedD.yn[j]);
			}
		}
		fflush(stdout); 
	
		switch (extra){
			/* extrapolation */
			case 0:	
				x_old[0] = forcedD.x[0];
				x_old[1] = forcedD.x[1];
				forcedD.gamma += (forcedD.dgamma / forcedD.extra);
				forcedD.k += (forcedD.dk / forcedD.extra);
				forcedD.B0 += (forcedD.dB0 / forcedD.extra);
				forcedD.B += (forcedD.dB / forcedD.extra);
				extra = 1;
				if (forcedD.extra == 0) extra = 2;
				break;
			case 1:
				forcedD.gamma += (forcedD.dgamma / forcedD.extra);
				forcedD.k += (forcedD.dk / forcedD.extra);
				forcedD.B0 += (forcedD.dB0 / forcedD.extra);
				forcedD.B += (forcedD.dB / forcedD.extra);
				if (k == extra - 1) extra = 2;
				break;
			case 2:
				{	double xtemp[2];
					xtemp[0] = forcedD.x[0];
					xtemp[1] = forcedD.x[1];
					forcedD.x[0] = 2.0 * xtemp[0] - x_old[0];
					forcedD.x[1] = 2.0 * xtemp[1] - x_old[1];
					x_old[0] = xtemp[0];
					x_old[1] = xtemp[1];
					forcedD.gamma += forcedD.dgamma;
					forcedD.k += forcedD.dk;
					forcedD.B0 += forcedD.dB0;
					forcedD.B += forcedD.dB;
				}
				break;
			}
	}
}


int fixed()
{
	int i,j,k;
	double prev[2];
	double delta;
	double dx, dy;

	i = 0;
	while (1){
		prev[0] = forcedD.x[0];
		prev[1] = forcedD.x[1];
		newton(forcedD);
		dx = forcedD.x[0] - prev[0];
		dy = forcedD.x[1] - prev[1];
		delta = (fabs(dx) + fabs(dy)) / 2.0;
#ifdef DEBUG
		printf("delta = %lf\n",delta);
#endif
		if (delta > forcedD.emax){
			fprintf(stderr, "\ndivergence.\n");
			exit(1);
		}
		i++;
		if (delta < forcedD.eps) break;
		if (i > forcedD.ite_max){
			printf("Sorry iteration > %d\n",forcedD.ite_max);
			exit(1);
		}
	}
	for (j = 0; j < 4; j++){ forcedD.z[j] = forcedD.x[j + 2]; }
	character(forcedD.z);
	return i;
}
	

int runge(n,h,x,t)
	int n;		/* number of the equations */
	double h;	/* strip */ 
	double x[];	/* variables */
	double t;	/* time */
{
	double k1[NN],k2[NN],k3[NN],k4[NN];
	double xtemp[NN];
	double hh;
	double thh;
	double th;
	int i;
	
	hh = h / 2.0;
	thh = t + hh;
	th = t + h;

	function(x,k1,t);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k1[i]*hh;
	function(xtemp,k2,thh);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k2[i]*hh;
	function(xtemp,k3,thh);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k3[i]*h;
	function(xtemp,k4,th);
	for (i=0; i<n; i++) x[i] += (k1[i]+2.0*(k2[i]+k3[i])+k4[i])*h/6.0;

}

int function(x,f,t)
	double x[],f[];
	double t;
{
	double p, q;

	p = -forcedD.gamma * exp(forcedD.gamma * x[0]);
	q = -forcedD.gamma * forcedD.gamma * exp(forcedD.gamma * x[0]);

	f[0] = -x[1] - exp(forcedD.gamma * x[0]);
	f[1] = x[0] + forcedD.k * x[1] - forcedD.B0 - forcedD.B * cos(t);

	f[2] = p * x[2] - x[3];				
	f[3] = x[2] + forcedD.k * x[3];
	f[4] = p * x[4] - x[5];				
	f[5] = x[4] + forcedD.k * x[5];
}



sysvar()
{
	int i,j,k;
	double t, h;

	h = M_PI*2.0/forcedD.m;

	forcedD.x[2] = 1.0; forcedD.x[3] = 0.0;
	forcedD.x[4] = 0.0; forcedD.x[5] = 1.0;

	for (i = 0; i < forcedD.l; i++){
		for (j = 0; j < forcedD.m; j++){
			t = h * (double)j;
			runge(NN,h,forcedD.x,t);
		}
		forcedD.xn[i] = forcedD.x[0];
		forcedD.yn[i] = forcedD.x[1];
	}
#ifdef DEBUG
    printf("a0[0][0] = %lf ", forcedD.x[2]);
	printf("a0[0][1] = %lf\n", forcedD.x[4]);
	printf("a0[1][0] = %lf ", forcedD.x[3]);
	printf("a0[1][1] = %lf\n", forcedD.x[5]);
#endif
}

int newton()
{
	int i;
	double prev[NN];
	double h[NN];
	double a[NN][NN];
	double b[NN];

	for ( i = 0; i < 2; i++) prev[i] = forcedD.x[i];

	sysvar();

	a[0][0] = forcedD.x[2] - 1.0;
	a[0][1] = forcedD.x[4];
	a[0][2] = prev[0] - forcedD.x[0]; 
	a[1][0] = forcedD.x[3];
	a[1][1] = forcedD.x[5] - 1.0;
	a[1][2] = prev[1] - forcedD.x[1]; 

	gauss(2,a,h);

	for ( i = 0; i < 2; i++) forcedD.x[i] = prev[i] + h[i];
}

int character(z)
	double z[];
{
	double alpha, beta, delta;
	double mu0, mu1, tr0, tr1;
	double rep;

	alpha = z[0] + z[3];
	beta = z[0] * z[3] - z[1] * z[2];
	delta = alpha * alpha - 4.0 * beta;

	if (delta > 0.0){
		mu0 = 0.5 * (alpha + sqrt(delta));
		mu1 = 0.5 * (alpha - sqrt(delta));
		tr0 = (mu0 - z[0])/z[2];
		tr1 = (mu1 - z[0])/z[2];
		z[0] = mu0;
		z[1] = mu1;
		z[2] = tr0;
		z[3] = tr1;
	} 
	else if (delta == 0.0){
		rep = 0.5 * alpha;
		z[0] = z[1] = rep;
		z[2] = z[3] = 0.0;
	} 
	else if (delta < 0.0) {
		z[0] = 0.5 * alpha; 			/* return real part */
		z[1] = 0.5 * sqrt(-1.0*delta); 	/* return imaginary part */
		z[2] = z[3] = 0.0;
	}
}

int gauss(n,a,h)
	int n;
	double a[NN][NN];
	double h[NN];
{
	int i,j,k;
	double sum;

	for (i = 0; i < n - 1; i++){
		pivoting(a, i, n);
		elimination(a, i, n);
	}

	for (i = n - 1; i >= 0; i--){
		sum = 0.0;
		for ( j = n -1; j > i; j--) sum += a[i][j]*h[j];
		h[i] = (a[i][n] - sum)/a[i][i];
	}
}


int pivoting(a,k,n)
	double a[NN][NN];
	int k, n;
{
	int i,j, posess = 0;
	double max = 0.0;
	double v0[NN],v1[NN];

	for (i = k; i < n; i++){
		double temp; 
		temp = fabs(a[i][k]); 	
		if (max < temp) { max = temp;	posess = i; }
	}
	for (i = 0; i < n + 1; i++){
		v0[i] = a[posess][i];
		v1[i] = a[k][i];
	}
	for (i = 0; i < n + 1; i++){
		a[posess][i] = v1[i];
		a[k][i] = v0[i];
	}
}

int elimination(a, k, n)
	double a[NN][NN];
	int k, n;
{
	double mp[NN][NN];
	int i,j;

	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (i == j) mp[i][j] = 1.0;
			else mp[i][j] = 0.0;
		}
	}
	for ( i = k + 1; i < n; i++) mp[i][k] = -a[i][k]/a[k][k];

	product(mp,a,n);
}

int product(mp,a,n)
	double mp[NN][NN];
	double a[NN][NN];
	int n;
{
	int i,j,k;
	double c[NN][NN];
	for (i = 0; i < n; i++){
		for (j = 0; j < n + 1; j++){
			c[i][j] = 0.0;
			for (k = 0; k < n; k++){
				c[i][j] += mp[i][k] * a[k][j];
			}
		}
	}
	for (i = 0; i < n; i++){
		for (j = 0; j < n+1; j++){
			a[i][j] = c[i][j];
		}
	}
}
