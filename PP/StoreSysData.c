#include "SysDefs.h"

SysData *StoreSysData(int argc, char **argv)
{
	int i, ii;
	FILE *fpin;
	SysData *sysdata;
    char cmd[BUFSIZ];

	if (argc == 1) { 
		fprintf(stderr,"Usage: %s datafile\n", argv[0]);
		return NULL;
	}

	if ((sysdata = (SysData *)calloc(sizeof (SysData), 1)) == NULL){
		fprintf(stderr, "cannot allocate memory\n");
		return NULL;
	}

    sprintf(cmd, "/usr/bin/cpp -P %s", argv[1]);

    if((fpin = popen(cmd, "r"))==NULL){
        fprintf(stderr,"cannot open %s\n", argv[1]);
        return NULL;
    }

	for (i = 0; i < NP; i++){
		fscanf(fpin,"%lf",&sysdata->param[i]);
	}

	fscanf(fpin,"%d",&sysdata->data);
	for (ii = 0; ii < sysdata->data; ii++){
		fscanf(fpin, "%d", &sysdata->line_wid[ii]);
		for (i = 0; i < NN; i++){ fscanf(fpin,"%lf",&sysdata->x0[ii][i]); }
		fscanf(fpin, "%lf", &sysdata->tsign[ii]);
		fscanf(fpin, "%lf", &sysdata->tend[ii]);
	}


	for (i = 0; i < NP; i++){
		fscanf(fpin,"%lf",&sysdata->dparam[i]);
	}

	fscanf(fpin,"%d",&sysdata->l);
	fscanf(fpin,"%d",&sysdata->m);
	fscanf(fpin,"%d",&sysdata->pflag);
	fscanf(fpin,"%d",&sysdata->write);

	xddp.line_wid = 1;
	xddp.line_att = SOLID;
	for (i=0;i<4;i++){ fscanf(fpin,"%d",&xddp.sv[i]); }
	for (i=0;i<4;i++){ fscanf(fpin,"%lf",&xddp.sc[i]); }
	for (i=0;i<2;i++){ fscanf(fpin,"%lf",&xddp.la[i]); }
	for (i=0;i<2;i++){ fscanf(fpin,"%d",&xddp.lav[i]); }
	for (i=0;i<2;i++){ fscanf(fpin,"%d",&xddp.fi[i]); }
	fscanf(fpin,"%d",&xddp.kse);
	fscanf(fpin,"%s",xddp.font);
	for (i=0;i<4;i++){ fscanf(fpin,"%d",&xddp.arrowtail[i]); }
	fscanf(fpin,"%d",&xddp.flag);

	if (sysdata->write){
		FILE *fp;

		if ((sysdata->write_fp = fopen(argv[2],"w"))==NULL){
				fprintf(stderr, "%s cannot open.\n",argv[2]);
				return NULL;
		}

		fp = sysdata->write_fp;
		fprintf(fp, "set terminal postscript  portrait plus \"Helvetica\" 18\n");
		fprintf(fp, "set size 1,0.7\n");
		fprintf(fp, "set ticslevel 0\n");
		fprintf(fp, "set output \"%s.ps\" \n", argv[2]);
		fprintf(fp, "set xlabel '\\size=24$x$'\n");
		fprintf(fp, "set ylabel '\\size=24$y$'\n");
		fprintf(fp, "set zlabel '\\size=24$z$'\n");
		fprintf(fp, "splot \"-\" notitle w l lt 1 lw 0.5\n");
		fflush(sysdata->write_fp);
	}
	xddparc.width = 4;
	xddparc.height = 4;
	xddparc.angle1 = 23040;
	xddparc.angle2 = 23040;
	xddparc.icol = 4;

	return sysdata;
}
