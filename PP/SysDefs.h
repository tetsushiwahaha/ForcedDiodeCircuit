#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <xhplot.h>

#define NN 2	/* dimension */
#define NUM 20
#define NP 4 	/* number of parameters */

#define GAMMA	0
#define K		1
#define B0		2
#define B		3

typedef struct sys_data { 
	double param[NP];
	double dparam[NP];

	int data;
	int line_wid[NUM];

	double h;
	int m;
	int l;

	double x0[NUM][NN];
	double tsign[NUM];
	double tend[NUM];

	int pflag;
	int write, graph;

	FILE *write_fp;

} SysData;

void DoEvent(double [], double, SysData *); 
double CalcSec(double [], double, double,  double, SysData *);
SysData *StoreSysData(int, char **);
void runge(int, double, double [], double, SysData *);

