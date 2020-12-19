#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define VE		18
#define PMAX	16

#define GAMMA	0
#define K		1
#define B0		2
#define B		3

#define PARAM	4 

struct FDstruct { 
	double param[PARAM];	/* prameters */
	double dparam[PARAM];	/* increments of parameters */
	
	double eps;			/* eps for newton's method */
	double emax;		/* error maximum boundary */

	int m;				/* integral strip */
	int l;				/* order of period */
	int extra;			/* extrapolation */
	int ite_max;		/* iteration maximum */
	int variable;		/* initial parameter number for VE */
	int increment;		/* initial incremental parameter */

	double x0[VE];		/* initial value */
	double x[VE];		/* state variables */
	double xn[PMAX];		/* periodic solution */
	double yn[PMAX];		/* periodic solution */
	double z[4];		/* characteristic variable */
	double mu;			/* over take eigen value */
	double grad;		/* gradient of shooting */
	int mode;
};

int main(int argc, char **argv)
{
	struct FDstruct *FD;
	int i, j, k = 0;
	int iteration, fcount;
	FILE *fopen(), *fpin, *fpout0;
	char name[BUFSIZ];
	char modestr[BUFSIZ];
	int fixed(struct FDstruct *);
	void nextstep(struct FDstruct *);

	if ((FD = calloc(1, sizeof(struct FDstruct))) == NULL){
		fprintf(stderr, "cannot allocate memory\n");
		exit(-1);
	}

	if (argc == 1){
		fprintf(stderr, "usage: %s datafile \n",argv[0]);
		exit(0);
	}

	sprintf(name, "/usr/bin/cpp -P %s", argv[1]);

	if((fpin = popen(name, "r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(-1);
	}

	sprintf(name, "%s.out",argv[1]); 
	if ((fpout0 = fopen(name,"w")) == NULL){
		fprintf(stderr, "cannot open %s\n",name);
		exit(-1);
	}

	/* read parameters for equation */

	for (i = 0; i < PARAM; i++){ fscanf(fpin,"%lf", &FD->param[i]); }
	for (i = 0; i < 2; i++){ fscanf(fpin,"%lf", &FD->x0[i]); }
	for (i = 0; i < PARAM; i++){ fscanf(fpin,"%lf", &FD->dparam[i]); }
	fscanf(fpin,"%s", modestr);
	fscanf(fpin,"%lf", &FD->eps);
	fscanf(fpin,"%lf", &FD->emax);
	fscanf(fpin,"%d", &FD->l);
	fscanf(fpin,"%d", &FD->m);
	fscanf(fpin,"%d", &FD->ite_max);
	fscanf(fpin,"%d", &FD->extra);
	fscanf(fpin,"%lf", &FD->grad);
	fscanf(fpin,"%d", &FD->increment);
	fscanf(fpin,"%d", &FD->variable);

	for (i = 0; i < 2; i++){ FD->x[i] = FD->x0[i]; }

	printf("i gamma  k    B0     B      x(1)     x(2)    mu1     mu2\n"); 

	while (1) {
		double mu;
		char s[BUFSIZ], s0[BUFSIZ];
		int compflag;
		int flg;
		int mode;

		if (!strcmp(modestr, "NS")) { FD->mode = 2; }
		else if (!strcmp(modestr, "PD") || !strcmp(modestr, "I")){ 
			FD->mode = 1;  mu = -1.0;
		} else if (!strcmp(modestr, "G") || !strcmp(modestr, "T")){ 
			FD->mode = 0;  mu = 1.0;
		}
		FD->mu = mu;

		iteration =	fixed(FD);

		if (fabs(FD->z[2]) < FD->eps){ 
			compflag = 1;
		} else { 
			compflag = 0;
		}
		sprintf(s,"%d %.5f %.5f %.5f %.5f %.7f %.7f", 
			iteration, FD->param[GAMMA], FD->param[K],
			FD->param[B0], FD->param[B], FD->x[0], FD->x[1]); 
		if (compflag){
			if (fabs(FD->x[2]+FD->x[5]) > 2.0){ flg = 'x'; } else { flg = 'o'; }
			sprintf(s0," %.3f+j%.3f (%.4f) %1c",FD->z[0], FD->z[1],
				sqrt(FD->z[0]*FD->z[0] + FD->z[1]*FD->z[1]), flg);
		} else {
			sprintf(s0," %.5f %.5f R", FD->z[0], FD->z[1]);
		}
		strncat(s, s0, BUFSIZ);
		if (FD->mode ==2 && !compflag) {
			printf("(false compulation...)\n");
		} else {
			puts(s);
			fputs(s, fpout0);
			fputs("\n", fpout0);
		}

		if ( FD->l != 1){
			for (j = 0; j < FD->l; j++){
				printf("     %.8f %.8f\n", FD->xn[j], FD->yn[j]);
				fprintf(fpout0, "     %.8f %.8f\n", 
					FD->xn[j], FD->yn[j]);
			}
			fflush(fpout0); 
		}
		fflush(stdout); 
		nextstep(FD);
	}
}


int fixed(FD)
	struct FDstruct *FD;
{
	int i,j,k;
	double prev[3];
	double delta;
	double dx, dy, dz;
	void newton(struct FDstruct *);
	void character(double *);

	i = 0;
	while (1){
		prev[0] = FD->x[0];
		prev[1] = FD->x[1];
		prev[2] = FD->param[FD->variable];
		newton(FD);
		dx = FD->x[0] - prev[0];
		dy = FD->x[1] - prev[1];
		dz = (FD->param[FD->variable] - prev[2])
			/ FD->param[FD->variable];
		if (fabs(FD->param[FD->variable]) < 0.01){
			dz = (FD->param[FD->variable] - prev[2]) / 0.01;
		}
		delta = (fabs(dx) + fabs(dy) + fabs(dz)) / 3.0;
#ifdef DEBUG
		printf("delta = %lf\n", delta); fflush(stdout);
#endif
		if (delta > FD->emax){
			fprintf(stderr, "Divergence.\n");
			exit(-1);
		}
		i++;
		if (delta < FD->eps) break;
		if (i > FD->ite_max){
			printf("Sorry iteration is over > %d\n",FD->ite_max);
			fflush(stdout);
			exit(-2);
		}
	}
	for (j = 0; j < 4; j++){ FD->z[j] = FD->x[j + 2]; }
	character(FD->z);
	return i;
}
	

void function(double *x, double *f, double t, struct FDstruct *FD)
{
	double p, q, gamma, k, b0, b;

	gamma = FD->param[GAMMA];
	k = FD->param[K];
	b0 = FD->param[B0];
	b = FD->param[B];

	p = -gamma * exp(gamma*x[0]);
	q = -gamma * gamma * exp(gamma * x[0]);

	f[ 0] = -x[1] - exp(gamma * x[0]);
	f[ 1] = x[0] + k * x[1] - b0 - b * cos(t);
	f[ 2] = p * x[ 2] - x[ 3];
	f[ 3] = x[ 2] + k * x[ 3];
	f[ 4] = p * x[ 4] - x[ 5];
	f[ 5] = x[ 4] + k * x[ 5];
	f[ 6] = p * x[ 6] - x[ 7];
	f[ 7] = x[ 6] + k * x[ 7];
	f[ 8] = p * x[ 8] - x[ 9] + q * x[2] * x[2];
	f[ 9] = x[ 8] + k * x[ 9] ;
	f[10] = p * x[10] - x[11] + q * x[2] * x[4];
	f[11] = x[10] + k * x[11];
	f[12] = p * x[12] - x[13] + q * x[4] * x[4];
	f[13] = x[12] + k * x[13];
	f[14] = p * x[14] - x[15] + q * x[2] * x[6];
	f[15] = x[14] + k * x[15];
	f[16] = p * x[16] - x[17] + q * x[4] * x[6];
	f[17] = x[16] + k * x[17];


	switch (FD->variable){
		case 0: /* gamma */
			f[6]  -= x[0] * exp(gamma * x[0]);
			f[14] -= (exp(gamma*x[0]) + gamma*x[0] * exp(gamma*x[0]))*x[2]; 
			f[16] -= (exp(gamma*x[0]) + gamma*x[0] * exp(gamma*x[0]))*x[4]; 
			break;
		case 1: /* k */
			f[7]  += x[1];
			f[15] += x[3];
			f[17] += x[5];
			break;
		case 2:
			f[7] -= 1.0;
			break;
		case 3:
			f[7] -= cos(t);
			break;
	}
}

void newton( struct FDstruct *FD)
{
	int i, j;
	double prev[VE];
	double h[VE];
	double a[VE][VE];
	double b[VE];
	double *x;
	int mode;
	double mu;

	void gauss(int, double [][VE], double[VE]);
	void sysvar(struct FDstruct *);

	mode = FD->mode;
	mu = FD->mu;

	x = FD->x;
	for ( i = 0; i < 2; i++) prev[i] = x[i];
	prev[2] = FD->param[FD->variable];
	sysvar(FD);

	/* make A = F'(x) | -F(x) to Newton Method for F'(x)dx = -F(x) */

	a[0][0] = x[2] - 1.0;
	a[0][1] = x[4];
	a[0][2] = x[6];
	a[0][3] = prev[0] - x[0]; 
	a[1][0] = x[3];
	a[1][1] = x[5] - 1.0;
	a[1][2] = x[7]; 
	a[1][3] = prev[1] - x[1]; 


#ifdef DEBUG
	printf("mode = %d\n", FD->mode);
#endif

	switch(FD->mode){
		case 0:
		case 1: 
			a[2][0] = -mu * (x[8] + x[11]) 
				+ x[ 8] * x[5] + x[2] * x[11] - x[10] * x[3] - x[4] * x[ 9];
			a[2][1] = -mu * (x[10] + x[13]) 
				+ x[10] * x[5] + x[2] * x[13] - x[12] * x[3] - x[4] * x[11];
			a[2][2] = -mu * (x[14] + x[17]) 
				+x[14] * x[5] + x[2] * x[17] - x[16] * x[3] - x[4] * x[15];
			a[2][3] = -(mu * mu - mu * (x[2] + x[5]) 
				+ x[2] * x[5] - x[4] * x[3]);
		break;
		case 2:
			a[2][0] = x[ 8] * x[5] + x[2] * x[11] - x[10] * x[3] - x[4] * x[ 9];
			a[2][1] = x[10] * x[5] + x[2] * x[13] - x[12] * x[3] - x[4] * x[11];
			a[2][2] = x[14] * x[5] + x[2] * x[17] - x[16] * x[3] - x[4] * x[15];
			a[2][3] = -(x[2] * x[5] - x[4] * x[3] - 1.0);
		break;
	}

#ifdef DEBUG
	for (i = 0; i < 3; i++){
		for (j = 0; j < 4; j++){
			printf("%.15f ", a[i][j]);
		}
		printf("\n");
	}
#endif

	gauss(3, a, h);

#ifdef DEBUG
	printf(">>>> %.15f %.15f %.15f\n", h[0], h[1], h[2]);
#endif

	for ( i = 0; i < 2; i++) x[i] = prev[i] + h[i];
	FD->param[FD->variable] = prev[2] + h[2];
}

void sysvar( struct FDstruct *FD)
{
	int i,j,k;
	double t, h;
	void runge(int, double, double *, double, struct FDstruct *);

	h = M_PI * 2.0 / FD->m;

	for (i = 2; i < VE; i++) FD->x[i] = 0.0;

	FD->x[2] = 1.0;  FD->x[5] = 1.0;

	for (i = 0; i < FD->l; i++){
		for (j = 0; j < FD->m; j++){
			t = h * (double)j;
			runge(VE, h, FD->x, t, FD);
		}
		FD->xn[i] = FD->x[0];
		FD->yn[i] = FD->x[1];
	}
}


void nextstep(struct FDstruct *d)
{
	static int count = 0, extra = 0; 
	int i, inc; static double x_old[3];
	double xx[2], absdp[PARAM], ptemp, grad, p, eps = 1.0e-5;
	void autopara(double *, struct FDstruct *, double *, int);

	switch (extra){ /* case 0 and 1 are extrapolation */
		case 0:		
			x_old[0] = d->x[0];
			x_old[1] = d->x[1];
			x_old[2] = d->param[d->variable];
			d->param[d->increment] += 
				(d->dparam[d->increment] / d->extra);
			extra++; count++;
			if (d->extra == 1) extra = 2;
			break;
		case 1:
			d->param[d->increment] 
				+= (d->dparam[d->increment] / d->extra);
			if (count == d->extra - 1){ extra = 2;}
			count++; 
			break;
		case 2:
			for (i = 0; i < PARAM; i++) {
				absdp[i] = fabs(d->dparam[i]);
			}
			grad = fabs(x_old[2] - d->param[d->variable]) 
				- d->grad * absdp[d->variable];
			if (grad > 0.0) {
				printf("change parameter \n");
				autopara(x_old, d, absdp, d->increment);
				extra = 0;
			} else {
				xx[0] = d->x[0]; xx[1] = d->x[1];
				ptemp = d->param[d->variable];
				d->x[0] = 2.0 * xx[0] - x_old[0];
				d->x[1] = 2.0 * xx[1] - x_old[1];
				d->param[d->variable] = 
					2.0 * (d->param[d->variable]) - x_old[2];
				x_old[0] = xx[0];
				x_old[1] = xx[1];
				x_old[2] = ptemp;
				d->param[d->increment] += d->dparam[d->increment];
				if (d->param[d->increment] < -eps) p = -0.5;
				else if (fabs(d->param[d->increment]) < eps) p = 0.0;
				else p = 0.5;
				d->param[d->increment] =  
					(int)(d->param[d->increment] / 
						absdp[d->increment] + p) * absdp[d->increment];
			}
			count = 0;
			break;
	}
}

void autopara(double *x_old, struct FDstruct *FD, 
	double *absdp, int inc)
{
	double av, p, eps = 1.0e-5; 

	FD->dparam[FD->variable] = absdp[FD->variable];
	if (x_old[2] > FD->param[FD->variable]){ 
		FD->dparam[FD->variable] = -FD->dparam[FD->variable]; }
	av = (x_old[2] + FD->param[FD->variable]) / 2.0;
	if (av < -eps) p = -0.5; 
	else if (fabs(av) < eps) p = 0.0;
	else p = 0.5;
	FD->param[FD->variable] = 
		(int)(av / absdp[FD->variable] + p) * absdp[FD->variable];
	printf("change to %lf %lf\n",
		 FD->param[FD->variable], FD->dparam[FD->variable]);
	fflush(stdout);
	FD->increment = FD->variable;	/* change parameter */
	FD->variable = inc;	/* change parameter */
}

/*********** DO NOT CHANGE BELOW **********************/

void character(double z[])
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

void gauss(int n, double a[VE][VE], double h[VE])
{
	int i,j,k;
	double sum;
	void pivoting(double [][VE], int, int);
	void elimination(double [][VE], int, int);

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

void pivoting(double a[VE][VE], int k, int n)
{
	int i,j, posess = 0;
	double max = 0.0;
	double v0[VE],v1[VE];

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

void elimination( double a[VE][VE], int k, int n)
{
	double mp[VE][VE];
	int i,j;
	void product(double [][VE], double [][VE], int);

	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (i == j) mp[i][j] = 1.0;
			else mp[i][j] = 0.0;
		}
	}
	for ( i = k + 1; i < n; i++) mp[i][k] = -a[i][k]/a[k][k];

	product(mp,a,n);
}

void product(double mp[VE][VE], double a[VE][VE], int n)
{
	int i,j,k;
	double c[VE][VE];
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

void runge(n,h,x,t,FD)
	int n;		/* number of the equations */
	double h;	/* strip */ 
	double x[];	/* variables */
	double t;	/* time */
	struct FDstruct *FD;
{
	double k1[VE],k2[VE],k3[VE],k4[VE];
	double xtemp[VE];
	double hh;
	double thh;
	double th;
	int i;
	void function(double *, double *, double, struct FDstruct *);
	
	hh = h / 2.0;
	thh = t + hh;
	th = t + h;

	function(x,k1,t,FD);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k1[i]*hh;
	function(xtemp,k2,thh,FD);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k2[i]*hh;
	function(xtemp,k3,thh,FD);
	for (i=0; i<n; i++) xtemp[i]=x[i]+k3[i]*h;
	function(xtemp,k4,th,FD);
	for (i=0; i<n; i++) x[i] += (k1[i]+2.0*(k2[i]+k3[i])+k4[i])*h/6.0;
}
