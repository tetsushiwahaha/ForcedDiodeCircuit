#include <X11/Xlib.h> 
#include <X11/Xutil.h>
#include <X11/keysym.h>
/* #include <floatingpoint.h> */
#include <stdio.h>
#include <math.h>
#ifdef BSD
#include <strings.h>
#endif BSD
#include <string.h>
#include "xhplot.h"
#define	DISPLAY	NULL
#define	NCOLORS	(16)
#define	MAXNPOINTS (4096)
#define	ARROWLONG (100)
#define	CHARLEN (255)
#define	STRINGLEN (255)
#define	XSIZE (1000)
#define	YSIZE (1000)
#define FALSE (0)
#define TRUE (! FALSE)
#define TICX (6)
#define TICY (6)
#define CFACT (2.5)
#define COFFX (200)
static Display *display;
static Window win;
static GC gc;
static XSetWindowAttributes winattr;
static unsigned long fg,bg;
static XColor color[NCOLORS];
static XFontStruct *font;
static char fontname[CHARLEN] = "vr-40";
static XPoint points[MAXNPOINTS];
static double xmag = 1.0,ymag = 1.0;
static int linestyle;
static double linelong;
static int xsize = XSIZE,ysize = YSIZE;
static reverse = FALSE;
static double logical_off=0.33;
#include "white.h"
#include "black.h"
#include "dashes.h"
#include "colors.h"
static XPoint arrowpoint[][6] = {
	{0,0}, {ARROWLONG,0}, {-2,5}, {23,-5}, {-23,-5}, {2,5},
	{0,0}, {0,-ARROWLONG}, {5,2}, {-5,-23}, {-5,23}, {5,-2}
};

int xhplot(mode,x,y,icol)
	int mode;
	double x,y;
	int icol;
{
	static int 
		sv[4] = { 130,830,50,750 },
		full[4] = { 0,XSIZE,0,YSIZE },
		lav[2] = { 5,5 },
		fi[2] = { 1,1 },
		kse = AXES,
		ix0,
		iy0,
		arrowtail[4] = { 420,840,20,460 },
		npoints,
		arcangle1 = 0,
		arcangle2 = 360*64,
		flag;
	static unsigned int 
		arcwidth = 8,
		archeight = 8;
	static double 
		xfactor,
		yfactor,
		sc[4] = { -2.0,2.0,-2.0,2.0 },
		la[2];
	static int ix,iy;
	static int 
		arrowver = HORIZONTAL,
		arrowhead = TAIL,
		arrowinv = NORMAL;
	static unsigned int xbmwidth,xbmheight;
	static Pixmap xbm;
	static int xbmxhot,xbmyhot;
	static char string[STRINGLEN] = "You must set a string.";
	switch (mode) {
	case WINDOW:
		full[1] = xsize = (int)x;
		full[3] = ysize = (int)y;
		return(0);
	case INIT_DRAW_RV:
		reverse = TRUE;
	case INIT_DRAW:
		init();
	case REDRAW:
		/* create default graphics context */
		gc = XCreateGC(display,win,NULL,NULL);
		XSetForeground(display,gc,fg);
		XSetBackground(display,gc,bg);
		data(sv,sc,&xfactor,&yfactor,la,lav,arrowtail);
		axes_data(la,lav,fi,&kse,arrowtail,&flag);
		axes(sv,lav,fi,&kse,&xfactor,&yfactor,sc,la,icol,flag);
		if (flag != NOTDRAW) {
			draw_arrow(HORIZONTAL,isgn(la[0]),arrowtail[0],arrowtail[1]);
			draw_arrow(VERTICAL,isgn(la[1]),arrowtail[2],arrowtail[3]);
		}
		cliprectangle(sv);
		break;
	
	case PSET:
		coordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		setcolor(icol);
		/* scalein(&ix,&iy); */
		putpixel(ix,iy);
		ix0 = ix; 
		iy0 = iy; 
		linelong = 0.0;
		break;
	case CYLPSET:
		fflush(stdout);
		cylcoordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		setcolor(icol);
		putpixel(ix,iy);
		ix0 = ix; 
		iy0 = iy; 
		linelong = 0.0;
		break;
	case PRELINE:
		coordinate(&ix0,&iy0,xfactor,yfactor,x,y,sv,sc);
		/* scalein(&ix0,&iy0); */
		linelong = 0.0;
		break;
	case LINE: 
		setcolor(icol);
		coordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		/* scalein(&ix,&iy); */
		line(ix0,iy0,ix,iy);
		ix0 = ix;
		iy0 = iy;
		break;
	case CYLLINE: 
		setcolor(icol);
		cylcoordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		line(ix0,iy0,ix,iy);
		ix0 = ix;
		iy0 = iy;
		break;
	case CLS: 
		XClearWindow(display,win);
		break;
	case LINEATT:
		setlineatt();
		break;
	case ARCSTRUCT:
		setarcstruct(&arcwidth,&archeight,&arcangle1,&arcangle2);
		break;
	case ARC:
		arcwidth = xddparc.width;
		archeight = xddparc.height;
		arcangle1 = xddparc.angle1;
		arcangle2 = xddparc.angle2;
		coordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		icol = xddparc.icol;
		setcolor(icol);
		arc(ix,iy,arcwidth,archeight,arcangle1,arcangle2);
		ix0 = ix; 
		iy0 = iy; 
		break;
	case INIT_SETP:
		npoints = 0;
		break;
	case SETPOINTS:
		coordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		scalein(&ix,&iy);
		setpoints(ix,iy,npoints++);
		break;
	case PSETS:
		setcolor(icol);
		putpixels(npoints);
		break;
	case LINES:
		setcolor(icol);
		lines(npoints);
		break;
	case FILL:
		fill(mode,npoints,icol);
		break;
	case FILLOpaque:
		fill(mode,npoints,icol);
		break;
	case CLIPMASK:
		clipmask(npoints);
		break;
	case CLIPFULL:
		cliprectangle(full);
		break;
	case FRAME:
		cliprectangle(full);
		flag = 1;
		struct2data(xddp,sv,sc,la,lav,fi,
			&kse,fontname,arrowtail,&flag);
		xfactor = ((double)(sv[1]-sv[0]))/(sc[1]-sc[0]);
		yfactor = ((double)(sv[3]-sv[2]))/(sc[3]-sc[2]);
		axes(sv,lav,fi,&kse,&xfactor,&yfactor,sc,la,icol,flag);
		if (flag != NOTDRAW) {
			draw_arrow(HORIZONTAL,isgn(la[0]),
				arrowtail[0],arrowtail[1]);
			draw_arrow(VERTICAL,isgn(la[1]),
				arrowtail[2],arrowtail[3]);
		}
		cliprectangle(sv);
		break;
	case CYLFRAME:
		cliprectangle(full);
		flag = 1;
		cylaxes(sv,lav,fi,&kse,&xfactor,&yfactor,sc,la,icol,flag);
		cliprectangle(sv);
		break;
	case DIRECT_RV:
		reverse = TRUE;
	case DIRECT:
		/* create default graphics context */
		init();
		linestyle=xddp.line_att;
		gc = XCreateGC(display,win,NULL,NULL);
		XSetForeground(display,gc,fg);
		XSetBackground(display,gc,bg);
		struct2data(xddp,sv,sc,la,lav,fi,
			&kse,fontname,arrowtail,&flag);		
		xfactor = ((double)(sv[1]-sv[0]))/(sc[1]-sc[0]);
		yfactor = ((double)(sv[3]-sv[2]))/(sc[3]-sc[2]);
		axes(sv,lav,fi,&kse,&xfactor,&yfactor,sc,la,icol,flag);
		if (flag != NOTDRAW) {
			draw_arrow(HORIZONTAL,isgn(la[0]),
				arrowtail[0],arrowtail[1]);
			draw_arrow(VERTICAL,isgn(la[1]),
				arrowtail[2],arrowtail[3]);
		}
		cliprectangle(sv);
		break;
	case CYLINDER:
		cliprectangle(full);
		cylaxes(sv,lav,fi,&kse,&xfactor,&yfactor,sc,la,icol,flag);
		cliprectangle(sv);
		break;
	case BITMAP:
		cliprectangle(full);
		setcolor(icol);
		putbitmap();
		break;
	case SETARROW:
		set_arrow(&arrowver,&arrowhead,&arrowinv);
		break;
	case ARROW:
		coordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		setcolor(icol);
		arrow(ix,iy,arrowver,arrowhead,arrowinv);
		break;
	case MKARROW:
		mkarrow();
		break;
	case SETBITMAP:
		readbitmap(&xbmwidth,&xbmheight,&xbm,&xbmxhot,&xbmyhot);
		break;
	case BITMAPS:
		coordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		setcolor(icol);
		scalein(&ix,&iy);
		bitmaps(ix,iy,xbmwidth,xbmheight,xbm,xbmxhot,xbmyhot);
		break;
	case SETFONT:
		setfont();
		break;
	case STRING:
		setcolor(icol);
		drawstring();
		break;
	case SETSTRING:
		setstring(string);
		break;
	case STRINGS:
		coordinate(&ix,&iy,xfactor,yfactor,x,y,sv,sc);
		setcolor(icol);
		scalein(&ix,&iy);
		drawstrings(ix,iy,string);
		break;
	case POINTER:
		pointer(xfactor,yfactor,sv,sc, &x, &y);
		break;
	case POINT:
		point(xfactor,yfactor,sv,sc, &x, &y);
		break;
	case CYLPOINT:
		cylpoint(xfactor,yfactor,sv,sc, &x, &y);
		break;
	case KEYIN:
		no_wait_key();
		break;
	case STATICKEYIN:
		key_in();
		break;
	case CLOSE:
		XCloseDisplay(display);
		return(0);
	case WAIT:
		waitseconds();
		return(0);
	case MAGNIFY:
		xmag = x;
		ymag = y;
		return(0);
	case BITMAPOUT:
		bitmapout();
		return(0);
	default:
		fprintf(stderr,"xhplot: %d is not a command.\n",mode);
		return(1);
	}
	XFlush(display);
	return(-1);
}
scalein(x,y)
int *x,*y;
{
	if (*x < 0) {
		*x = 0;
		fprintf(stderr,"<");
	} 
	else if (*x > (int)((double)xsize/xmag)) {
		*x = (int)((double)xsize/xmag);
		fprintf(stderr,">");
	}
	if (*y < 0) {
		*y = 0;
		fprintf(stderr,"^");
	}
	else if (*y > (int)((double)ysize/ymag)) {
		*y = (int)((double)ysize/ymag);
		fprintf(stderr,"_");
	}
}
waitseconds()
{
	unsigned seconds;
#ifdef VERBOSE
	printf("sleeping seconds\n");
#endif VERBOSE
	scanf("%d",&seconds);
	sleep(seconds);
}

pointer(xfactor,yfactor,sv,sc,x,y)
	int sv[];
	double xfactor,yfactor,sc[];
	double *x,*y;
{
	int ix,iy;
	int button;
	XEvent event;
	XSelectInput(display,win,ButtonPressMask);
	for (;;) {
		XNextEvent(display,&event);
		if (event.type == ButtonPress) {
			if ((button = event.xbutton.button) == 3) break;
			ix = event.xbutton.x;
			iy = event.xbutton.y;
			switch (button) {
			case 1: 
				fprintf(stderr,"%d %d\n",ix,iy);
				break;
			case 2:
				arccoordinate(ix,iy,xfactor,yfactor,x,y,sv,sc);
				fprintf(stdout,"%f %f\n",*x,*y);
				fflush(stdout);
				break;
			default:
				fprintf(stderr,"%d: no such button supported\n",
					button);
			}
		}
	}
	xddpret.x = *x;
	xddpret.y = *y;
}

point(xfactor,yfactor,sv,sc,x,y)
	int sv[];
	double xfactor,yfactor,sc[];
	double *x,*y;
{
	int ix,iy;
	int button;
	Bool test;
	XEvent event;
	XSelectInput(display,win,ButtonPressMask|KeyPressMask);
	test = XCheckMaskEvent(display,ButtonPressMask,&event);
	if (test == True){
		button = event.xbutton.button;
		if (event.type == ButtonPress) {
			ix = event.xbutton.x;
			iy = event.xbutton.y;
			switch (button) {
			case 1: 
					arccoordinate(ix,iy,xfactor,yfactor,x,y,sv,sc);
					fprintf(stderr,"%f %f\n",*x,*y);
					break;
			}
		}
		xddpret.x = *x;
		xddpret.y = *y;
		xddpret.key = 0;
	} else {
		xddpret.key = -1;
	}
}

cylpoint(xfactor,yfactor,sv,sc,x,y)
	int sv[];
	double xfactor,yfactor,sc[];
	double *x,*y;
{
	int ix,iy;
	int button;
	Bool test;
	XEvent event;
	XSelectInput(display,win,ButtonPressMask|KeyPressMask);
	test = XCheckMaskEvent(display,ButtonPressMask,&event);
	if (test == True){
		button = event.xbutton.button;
		if (event.type == ButtonPress) {
			ix = event.xbutton.x;
			iy = event.xbutton.y;
			switch (button) {
			case 1: 
					arccylcoordinate(ix,iy,xfactor,yfactor,x,y,sv,sc);
					fprintf(stderr,"%f %f\n",*x,*y);
					break;
			}
		}
	}	
	xddpret.x = *x;
	xddpret.y = *y;
}

cliprectangle(cp)
int cp[];
{
	double magx(),magy();
	XRectangle rectangles;
	rectangles.x = (short)magx(cp[0]);
	rectangles.y = (short)magy(cp[2]);
	rectangles.width  = (unsigned short)(magx(cp[1]-cp[0]));
	rectangles.height = (unsigned short)(magy(cp[3]-cp[2]));
	XSetClipRectangles(display,gc,0,0,&rectangles,1,Unsorted);
}

setpoints(ix,iy,n)
int ix,iy,n;
{
	double magx(),magy();
	if (n >= MAXNPOINTS) {
		fprintf(stderr,"%d is maximum data number\n",n);
		exit(! NULL);
	}
	points[n].x = magx(ix);
	points[n].y = magy(iy);
}

clipmask(n)
int n;
{
	Pixmap bclip;
	GC gclip;
	linkpoints(&n);
	bclip = XCreatePixmap(display,RootWindow(display,0),
		xsize,ysize,1);
	gclip = XCreateGC(display,bclip,0,0);
	XSetFillStyle(display,gclip,FillSolid); 
	XSetFillRule(display,gclip,WindingRule);
	XSetForeground(display,gclip,0L);
	XFillRectangle(display,bclip,gclip,0,0,(unsigned int)xsize,(unsigned int)ysize);
	XSetForeground(display,gclip,1L);
	XFillPolygon(display,bclip,gclip,
		points,n,Complex,CoordModeOrigin);
	XSetClipMask(display,gc,bclip);
	XFreePixmap(display,bclip);
}
linkpoints(n)
int *n;
{
	if ((points[*n-1].x !=  points[0].x)
	|| (points[*n-1].y !=  points[0].y)) {
		++ *n;
		points[*n-1].x = points[0].x;
		points[*n-1].y = points[0].y;
	}
}

fill(flag,n,icol)
int flag,n,icol;
{
	unsigned int width,height;
	Pixmap bitmap;
	int xhot,yhot;
	switch (flag) {
	case FILL: 
		XSetFillStyle(display,gc,FillStippled); 
		break;
	case FILLOpaque: 
		XSetFillStyle(display,gc,FillOpaqueStippled); 
		break;
	}
	linkpoints(&n);
	setcolor(icol);
	readbitmap(&width,&height,&bitmap,&xhot,&yhot);
	XSetStipple(display,gc,bitmap);
	XSetFillRule(display,gc,WindingRule);
	XFillPolygon(display,win,gc,
		points,n,Complex,CoordModeOrigin);
	if (reverse)
		XSetStipple(display,gc,XCreateBitmapFromData(display,win,white_bits,white_width,white_height));
	else
		XSetStipple(display,gc,XCreateBitmapFromData(display,win,black_bits,black_width,black_height));
}
putbitmap()
{
	double magx(),magy();
	unsigned int width,height;
	Pixmap bitmap;
	int xhot,yhot,xdest,ydest;
	readbitmap(&width,&height,&bitmap,&xhot,&yhot);
#ifdef VERBOSE
	printf("position: ");
#endif VERBOSE
	scanf("%d %d",&xdest,&ydest);
	XCopyPlane(display,bitmap,win,gc,0,0,width,height,(unsigned int)magx(xdest),(unsigned int)magy(ydest),1L);
}

bitmaps(xdest,ydest,width,height,bitmap,xhot,yhot)
unsigned int width,height;
Pixmap bitmap;
int xhot,yhot,xdest,ydest;
{
	double magx(),magy();
	xhot = yhot = 0;
	XCopyPlane(display,bitmap,win,gc,xhot,yhot,width,height,(unsigned int)(magx(xdest)-width/2),(unsigned int)(magy(ydest)-height/2),1L);
}

readbitmap(width,height,bitmap,xhot,yhot)
unsigned int *width,*height;
Pixmap *bitmap;
int *xhot,*yhot;
{
	char bitmapfile[CHARLEN];
	if (XReadBitmapFile(display,win,bitmapfile,
		width,height,bitmap,xhot,yhot) != BitmapSuccess) {
		fprintf(stderr,"%s? No such bitmap file.\n",bitmapfile);
		exit(! NULL);
	}
}

setlineatt()
{
#ifdef VERBOSE
	printf("line attributes\n ");
	printf("width style: ");
#endif VERBOSE
/*
	if (
	scanf("%d %d",&width,&style) 
	!= 2) getchar();
*/
	switch (xddp.line_att) {
	case SOLID:
		linestyle = xddp.line_att;
		XSetLineAttributes(display,gc,
			xddp.line_wid,LineSolid,CapRound,JoinRound);
		break;
	case DASH:
	case DOTTED:
	case DOT_DASHED:
	case SHORT_DASHED:
	case LONG_DASHED:
	case ODD_DASHED:
		linestyle = xddp.line_att;
		XSetLineAttributes(display,gc,xddp.line_wid,
			LineOnOffDash,CapButt,JoinRound);
		break;
	}
}

coordinate(ix,iy,xfactor,yfactor,x,y,sv,sc)
int *ix,*iy,sv[];
double xfactor,yfactor,x,y,sc[];
{

	*ix = (int)(xfactor*(x-sc[0]))+sv[0];
	*iy = sv[3]-(int)(yfactor*(y-sc[2]));
}

cylcoordinate(ix,iy,xfactor,yfactor,x,y,sv,sc)
	int *ix,*iy,sv[];
	double xfactor,yfactor,x,y,sc[];
{
	double delta0,delta1;
	double sigma;

	delta0 = (sv[1] - sv[0])/2.0;
	delta1 = (sv[3] - sv[2])/2.0;

	sigma = (y - sc[2]) / (sc[3] - sc[2]) * delta0 * (1.0 - logical_off)
		+ delta0 * logical_off;
	*ix = (int) (sigma * cos(x)) + delta0 + sv[0];
	*iy = sv[3] - delta1 - (int) (sigma * sin(x));

}

arccoordinate(ix,iy,xfactor,yfactor,x,y,sv,sc)
	int ix,iy,sv[];
	double xfactor,yfactor,*x,*y,sc[];
{
	double arcmagx(),arcmagy();
	*x = sc[0] + (arcmagx(ix) - (double)sv[0])/xfactor;
	*y = sc[2] - (arcmagy(iy) - (double)sv[3])/yfactor;
}

arccylcoordinate(ix,iy,xfactor,yfactor,x,y,sv,sc)
	int ix,iy,sv[];
	double xfactor,yfactor,*x,*y,sc[];
{
	double delta0,delta1;
	double sigma;


	*x = acos((ix - delta0 - sv[0])/sigma);
	*y = asin((-iy + sv[3] - delta1)/sigma);

}
putpixel(ix,iy)
int ix,iy;
{
	double magx(),magy();
	XDrawPoint(display,win,gc,(int)magx(ix),(int)magy(iy));
}
putpixels(n)
int n;
{
	XDrawPoints(display,win,gc,points,n,CoordModeOrigin);
}
line(ix0,iy0,ix,iy)
	int ix0,iy0,ix,iy;
{
	double magx(),magy();
	switch (linestyle) {
	case SOLID:
		break;
	case DASH:
	case DOTTED:
	case DOT_DASHED:
	case SHORT_DASHED:
	case LONG_DASHED:
	case ODD_DASHED:
		linelong += sqrt( 
			(magx(ix)-magx(ix0))*(magx(ix)-magx(ix0)) + 
			(magy(iy)-magy(iy0))*(magy(iy)-magy(iy0)) );
		XSetDashes(display,gc
			,((int)linelong) % dash_total_size[linestyle-DASH]
			,dash_list[linestyle-DASH]
			,dash_list_length[linestyle-DASH]);
		break;
	}
	XDrawLine(display,win,gc,(int)magx(ix0),(int)magy(iy0),(int)magx(ix),(int)magy(iy));
}

lines(n)
	int n;
{
	XDrawLines(display,win,gc,points,n,CoordModeOrigin);
}

int setarcstruct(width,height,angle1,angle2)
unsigned int *width,*height;
int *angle1,*angle2;
{
#ifdef VERBOSE
	printf("arc structure\n ");
	printf("width height angle1 angle2(degree*64): ");
#endif VERBOSE
/*
	if (
	scanf("%d %d %d %d",width,height,angle1,angle2) 
	!= 4) getchar();
	*/
	*width = xddparc.width;
	*height = xddparc.height;
	*angle1 = xddparc.angle1;
	*angle2 = xddparc.angle2;
}

arc(ix,iy,width,height,angle1,angle2)
int ix,iy,angle1,angle2;
unsigned int width,height;
{
	double magx(),magy();
	ix = magx((int)(ix - width/2));
	iy = magy((int)(iy - height/2));
	width = magx((int)width);
	height = magy((int)height);
	/* Bug:
	XDrawArc(display,win,gc,ix,iy,width,height,angle1,angle2);
	*/
	XSetArcMode(display,gc,ArcPieSlice);
	XFillArc(display,win,gc,ix,iy,width,height,angle1,angle2);
}

setcolor(icol)
int icol;
{
	/* enter color code into a gc */
	XSetForeground(display,gc,color[icol].pixel);
}

data(sv,sc,xfactor,yfactor,la,lav,arrowtail)
int sv[],lav[],arrowtail[];
double sc[],*xfactor,*yfactor,la[];
{
#ifdef VERBOSE
	printf("HK-Graph Ver.3.2 on X11R4\n");
	printf("view [%d %d %d %d]: ",sv[0],sv[1],sv[2],sv[3]);
#endif VERBOSE
	if (
	scanf("%d %d %d %d",&sv[0],&sv[1],&sv[2],&sv[3]) 
	!= 4) getchar();
#ifdef VERBOSE
	printf("scale [%f %f %f %f]: ",sc[0],sc[1],sc[2],sc[3]);
#endif VERBOSE
	if (
	scanf("%lf %lf %lf %lf",&sc[0],&sc[1],&sc[2],&sc[3]) 
	!= 4) getchar();
	set_misc(sv,sc,xfactor,yfactor,la,lav,arrowtail);
}

axes_data(la,lav,fi,kse,arrowtail,flag)
double la[];
int lav[],fi[],*kse,arrowtail[];
int *flag;
{
	do {
#ifdef VERBOSE
		printf("axes [%f %f %d %d]: ",la[0],la[1],lav[0],lav[1]);
#endif VERBOSE
		if (
		scanf("%lf %lf %d %d", &la[0],&la[1],&lav[0],&lav[1]) 
		!= 4) getchar();
		fi[0] = fiest(la[0]*(double)lav[0]);
		fi[1] = fiest(la[1]*(double)lav[1]);
#ifdef VERBOSE
		printf("fixed [%d %d]: ",fi[0],fi[1]);
#endif VERBOSE
		if (
		scanf("%d %d",&fi[0],&fi[1]) 
		!= 2) getchar();
#ifdef VERBOSE
		printf("AXES(1) or GRID(2) [%d]: ",*kse);
#endif VERBOSE
		if (
		scanf("%d",kse) 
		!= 1) getchar();
		setfont();
#ifdef VERBOSE
		printf("arrows [%d %d %d %d]: ",arrowtail[0],arrowtail[1],arrowtail[2],arrowtail[3]);
#endif VERBOSE
		if (
		scanf("%d %d %d %d",&arrowtail[0],&arrowtail[1],&arrowtail[2],&arrowtail[3]) 
		!= 4) getchar();
		*flag = 1;
#ifdef VERBOSE
		printf("OKDRAW(1), NOTDRAW(2) or SETAGAIN(3) [%d]: ",*flag);
#endif VERBOSE
		if (
		scanf("%d",flag)
		!= 1) getchar();
	} while (*flag == SETAGAIN);
}

init()
{
	char tmp[256];
	int i;
	XEvent event;
	XColor exact_color;
	XClassHint class_hints;
	if ((display = XOpenDisplay(DISPLAY)) == NULL) {
		fprintf(stderr,"unable to open display\n");
		exit(! NULL);
	}
	/* set foreground and background colors */
	if (reverse) {
		fg = BlackPixel(display,0);
		bg = WhitePixel(display,0);
		for (i = 0; i < 4; i++) {
			strcpy(tmp,colornames[i]);
			strcpy(colornames[i],colornames[7-i]);
			strcpy(colornames[7-i],tmp);
		}
	} else {
		fg = WhitePixel(display,0);
		bg = BlackPixel(display,0);
	}
	/* create a window */
	winattr.background_pixel = bg;
	winattr.border_pixel = fg;
	winattr.event_mask = ButtonPressMask|ExposureMask|StructureNotifyMask;
        winattr.bit_gravity = CenterGravity;
        winattr.backing_store = Always;
	win = XCreateWindow(display,RootWindow(display,0),
	0,0,(unsigned int)xsize,(unsigned int)ysize,0, /* location and size */
	CopyFromParent,CopyFromParent,CopyFromParent, /* depth,class,visual */
	CWBackPixel|CWBorderPixel|CWEventMask|CWBackingStore,&winattr);
	XStoreName(display,win,
	class_hints.res_name = "xhplot"	);
        class_hints.res_class = "XHplot";
        XSetClassHint(display,win,&class_hints);
	XMapWindow(display,win);
	/* wait until this window will be mapped */
	XMaskEvent(display,StructureNotifyMask,&event);
	/* get pixel values */
	for (i = 0; i < NCOLORS; i++)
		XAllocNamedColor(display,DefaultColormap(display,0),
		colornames[i], /* specify the color with string */
		&(color[i]),&exact_color); /* get color structs */
}

axes(sv,lav,fi,kse,xfactor,yfactor,sc,la,icol,flag)
int sv[],lav[],fi[],*kse,flag;
double *xfactor,*yfactor,sc[],la[];
int icol;
{
	static int ticx = TICX,ticy = TICY;
	int ix,iy,k,lt0,lt1,lt2,flag_ur,flag_c;
	double t;
	double sgn();
	setcolor(icol);
	switch (*kse) {
	case GRID:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = sv[3];
		lt1 = lt0-ticy;
		lt2 = sv[2]; 
		break;
	case AXES_UR:
	case AXES_U:
		flag_ur = TRUE;
		flag_c = FALSE;
		lt0 = sv[2];
		lt1 = lt0+ticy;
		lt2 = lt0+ticy*2;
		break;
	case AXES_C:
		flag_ur = FALSE;
		flag_c = TRUE;
		lt0 = (sv[3]+sv[2])/2;
		lt1 = lt0-ticy;
		lt2 = lt0-ticy*2;
		break;
	case ABNORMAL:
		ticy = -ticy;
		break;
	case FRAMEONLY:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = sv[3];
		lt1 = lt0;
		lt2 = lt0;
		break;
	default:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = sv[3];
		lt1 = lt0-ticy;
		lt2 = lt0-ticy*2;
		break;
	}
	t = sc[0]; 
	k = 0;
	while (sgn(sc[1]-sc[0])*(t - (sc[1]+0.01*la[0])) <= 0.0){
		ix = (int)(*xfactor*(t-sc[0]))+sv[0];
		if (k%lav[0] == 0){
			if (flag != NOTDRAW && flag != NOLINE) 
				line(ix,lt0,ix,lt2);
			if (flag != NOTDRAW && flag != NOMOJI) 
				moji(HORIZONTAL,fi[0],ix,lt0,t,flag_ur);
			k = 0;
		} else {
			if (flag != NOTDRAW && flag != NOLINE) 
				line(ix,lt0,ix,lt1);
		}
		t += la[0];
		k++;
	}
	if (flag_c) {
		if (flag != NOTDRAW && flag != NOLINE) 
			line(sv[0],lt0,sv[1],lt0);
	} else {
		if (flag != NOTDRAW && flag != NOLINE) {
			line(sv[0],sv[3],sv[1],sv[3]);
		}
	}
	switch (*kse) {
	case GRID:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = sv[0];
		lt1 = lt0+ticx;
		lt2 = sv[1]; 
		break;
	case AXES_UR:
	case AXES_R:
		flag_ur = TRUE;
		flag_c = FALSE;
		lt0 = sv[1];
		lt1 = lt0-ticx;
		lt2 = lt0-ticx*2;
		break;
	case AXES_C:
		flag_ur = FALSE;
		flag_c = TRUE;
		lt0 = (sv[1]+sv[0])/2;
		lt1 = lt0+ticy;
		lt2 = lt0+ticy*2;
		break;
	case ABNORMAL:
		ticx = -ticx;
		break;
	case FRAMEONLY:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = sv[0];
		lt1 = sv[0];
		lt2 = sv[0];
		break;
	default:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = sv[0];
		lt1 = lt0+ticx;
		lt2 = lt0+ticx*2;
		break;
	}
	t = sc[2]; 
	k = 0;
	while (sgn(sc[3]-sc[2])*(t - (sc[3]+0.01*la[1])) <= 0.0){
		iy = sv[3]-(int)(*yfactor*(t-sc[2]));
		if (k%lav[1] == 0) {
			if (flag != NOTDRAW && flag != NOLINE) 
				line(lt0,iy,lt2,iy);
			if (flag != NOTDRAW && flag != NOMOJI) 
				moji(VERTICAL,fi[1],iy,lt0,t,flag_ur);
			k = 0;
		} else {
			if (flag != NOTDRAW && flag != NOLINE) 
				line(lt0,iy,lt1,iy);
		}
		t += la[1];
		k++;
	}
	if (flag_c) {
		if (flag != NOTDRAW && flag != NOLINE) 
			line(lt0,sv[3],lt0,sv[2]);
	} else {
		if (flag != NOTDRAW && flag != NOLINE) {
			line(sv[0],sv[3],sv[0],sv[2]);
		}
	}
	if (flag != NOTDRAW && flag != NOLINE && flag_c != TRUE) {
		line(sv[0],sv[2],sv[1],sv[2]);
		line(sv[1],sv[3],sv[1],sv[2]);
	}
}

cylaxes(sv,lav,fi,kse,xfactor,yfactor,sc,la,icol,flag)
	int sv[],lav[],fi[],*kse,flag;
	double *xfactor,*yfactor,sc[],la[];
	int icol;
{
	int center_x, center_y, width, height;
	int start_x, start_y;
	double xfact;
	static int ticx = TICX, ticy = TICY;
	int ix,iy,k,lt0,lt1,lt2,flag_ur,flag_c;
	double t, strip;
	double sgn();
	double angle;

	setcolor(icol);

	width = (sv[1]-sv[0]);
	height = (sv[3]-sv[2]);
	start_x = sv[0] + (int)(width * 2.0 / 3.0);
	start_y = sv[2] + (int)(height / 2.0);
	xfact = ( width / (3.0 * (sc[3] - sc[2])));
	strip = 2*M_PI / lav[0];

	switch (*kse) {
	default:
		flag_ur = TRUE;
		flag_c = FALSE;
		lt0 = start_y;
		lt1 = lt0-ticy;
		lt2 = lt0-ticy*2;
		break;
	}
	for ( k = 0; k < lav[0]; k++ ){
		int off_x, off_y;
		double s,c;

		angle = strip * k;
		s = sin(angle);
		c = cos(angle);
		off_x = 2 * (int)(ticx * c);
		off_y = 2 * (int)(ticy * s);
		ix = (int)(sv[0] + width / 2.0 * (1 + c));
		iy = (int)(sv[2] + height /2.0 * (1 + s));
		line(ix, iy, ix + off_x, iy + off_y);
		/*
		moji(HORIZONTAL, fi[0],
			ix + off_x*2, iy + off_y * 2, angle, flag_c);
		*/
	}
	XDrawArc(display, win, gc, sv[0] + width / 3, sv[2] + height / 3,
		width / 3, height / 3, 0, 360*64);
	XDrawArc(display, win, gc, sv[0], sv[2], width, height, 0, 360*64);

	switch (*kse) {
	case GRID:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = sv[3];
		lt1 = lt0-ticy;
		lt2 = sv[2]; 
		break;
	case AXES_UR:
	case AXES_U:
		flag_ur = TRUE;
		flag_c = FALSE;
		lt0 = sv[2];
		lt1 = lt0+ticy;
		lt2 = lt0+ticy*2;
		break;
	case AXES_C:
		flag_ur = FALSE;
		flag_c = TRUE;
		lt0 = (sv[3]+sv[2])/2;
		lt1 = lt0-ticy;
		lt2 = lt0-ticy*2;
		break;
	case ABNORMAL:
		ticy = -ticy;
	default:
		flag_ur = FALSE;
		flag_c = FALSE;
		lt0 = start_y;
		lt1 = lt0-ticy;
		lt2 = lt0-ticy*2;
		break;
	}
#if DEBUG 
	printf("lt0 = %d, lt1 = %d, lt2 = %d\n",lt0,lt1,lt2);
	fflush(stdout);
#endif
	t = sc[2]; 
	k = 0;
	while (sgn(sc[3] - sc[2])*(t - sc[3]) <= 0.0){
		ix = (int)(xfact * (t - sc[2])) + start_x;
		if (k%lav[1] == 0 ) {
				line(ix,lt0,ix,lt2);
				moji(HORIZONTAL, fi[0], ix, lt0, t, flag_ur);
				k = 0;
		} else {
				line(ix,lt0,ix,lt1);
		}
		t += la[1];
		k++;
	}
	line(start_x, start_y, start_x + width/6*2  ,start_y);
}

double 
sgn(x)
double x;
{
	if (x >= 0.0)
		return( 1.0);
	else
		return(-1.0);
}
gfcvt(val,ndigit,s)
double val;
int ndigit;
char s[];
{
#ifdef SUN
	char str[CHARLEN];
	char format[CHARLEN];
	int i = 0,j = 0,
	int dec, ndec, sign, ndig;
	if (fabs(val) * pow(10.0,(double)ndigit) < .9) {
		s[0] = '0'; s[1] = '\0';
		return(-1);
	}

        if (sign == 1) s[i++] = '-';
	if (dec <= 0) s[i++] = '0';
        if (dec < 0) {
                s[i++] = '.';
                ndec = dec;
		ndig = ndigit;
                while (((ndec++) != 0)&&((ndig--) != 0)) s[i++] = '0';
        }
        do {
                if ((j == dec)&&(ndigit != 0)) s[i++] = '.';
	} while ((s[i++] = str[j++]) != '\0');
#else
	char format[BUFSIZ];
	sprintf(format, "%%.%df", ndigit);
	sprintf(s, format, val); 
#endif
}
moji(ise,jfx,lbl,lt,t,flag_ur)
int ise,jfx,lbl,lt;
double t;
{
	double magx(),magy();
	char string[31];
	int len,width,height;
	XCharStruct cs;
	int dir,ascent,descent,x,y,gap;

	gfcvt(t, jfx, string);
	len = strlen(string);
	XTextExtents(font,string,len,&dir,&ascent,&descent,&cs);
	width = cs.rbearing-cs.lbearing;
	height = cs.ascent-cs.descent;
	switch (ise) {
	case HORIZONTAL:
		x = magx(lbl)-(int)((double)width/2.0);
		y = magy(lt);
		gap = (int)((double)height*magx(10)/10.0);
		if (flag_ur)
			y -= gap;
		else
			y += (height + gap);
		break;
	case VERTICAL:
		x = magx(lt);
		gap = (int)((double)height*magx(6)/10.0);
		if (flag_ur)
			x += gap;
		else
			x -= (width + gap);
		y = magy(lbl)+(int)((double)height/2.0);
		break;
	}
	XDrawString(display,win,gc,x,y,string,len);
} 

set_arrow(ver,head,inv)
int *ver,*head,*inv;
{
	int arrowlong;
#ifdef VERBOSE
	printf("horizontal/vertical tail/head normal/inverse long: ");
#endif VERBOSE
	scanf("%d %d %d %d",ver,head,inv,&arrowlong);
	arrowpoint[0][1].x =  arrowlong;
	arrowpoint[1][1].y = -arrowlong;
}
arrow(ix,iy,ver,head,inv)
int ix,iy,ver,head,inv;
{
	int s;
	s = 1-2*inv;
	if (head == HEAD) {
		if (ver == HORIZONTAL)
			ix -= s * (arrowpoint[0][1].x + arrowpoint[0][2].x + arrowpoint[0][3].x);
		if (ver == VERTICAL)
			iy -= s * (arrowpoint[1][1].y + arrowpoint[1][2].y + arrowpoint[1][3].y);
	}
	draw_arrow(ver,s,ix,iy);
}
draw_arrow(i,s,ix,iy)
int i,s,ix,iy;
{
	double magx(),magy();
	XPoint p[2][6];
	int j;
	XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinMiter);
	XSetFillRule(display,gc,WindingRule);
	p[i][0].x = magx(ix);
	p[i][0].y = magy(iy);
	for (j = 1; j < 6; j++) {
		p[i][j].x = s*magx((int)arrowpoint[i][j].x);
		p[i][j].y = s*magy((int)arrowpoint[i][j].y);
	}
	XDrawLines(display,win,gc,p[i],6,CoordModePrevious);
	XFillPolygon(display,win,gc,p[i],6,Complex,CoordModePrevious);
}
double 
magx(x)
int x;
{
	return((double)x * xmag);
}
double 
magy(y)
int y;
{
	return((double)y * ymag);
}
double 
arcmagx(x)
int x;
{
	return((double)x / xmag);
}
double 
arcmagy(y)
int y;
{
	return((double)y / ymag);
}
int 
isgn(x)
double(x);
{
	if (x >= 0.0)
		return(1);
	else
		return(-1);
}
int
fiest(x)
double x;
{
	double pow(),xpow;
	int i;
	unsigned long n;
	for (i = 0; i < 8; i++)
		if ((n = (unsigned long)((xpow = x * pow(10.0,(double)i)) + 0.5) * 1000) != 0 && n == (unsigned long)(xpow * 1000.0 + 0.5))
			break;
	return (i);
}
setfont()
{
	char fontname0[CHARLEN];
	int i;
	i = 0;
	for (;;) {
#ifdef VERBOSE
		printf("font [%s]: ",fontname);
#endif VERBOSE
		scanf("%s",fontname0);
		if (strchr(fontname0,'/') == NULL)
			strcpy(fontname,fontname0);
		if ((font = XLoadQueryFont(display,fontname)) == NULL) {
			fprintf(stderr,"%s? No such font.\n",fontname);
			if (++i == 8) exit(1);
		} else
			break;
	}
	/* enter the font into the graphic context */
	XSetFont(display,gc,font->fid);
}
drawstring()
{
	int x,y;
	char string[STRINGLEN];
	double magx(),magy();
	setstring(string);
#ifdef VERBOSE
	printf("position [%d %d]: ",x,y);
#endif VERBOSE
	scanf("%d %d",&x,&y);
	XDrawString(display,win,gc,(int)magx(x),(int)magy(y),string,(int)strlen(string));
}
setstring(string)
char string[];
{
	int i,s;
#ifdef VERBOSE
	printf("string [%s]: ",string);
#endif VERBOSE
	scanf("%s",string);
	if (string[0] == '$' && string[1] != '$') {
		i = 0;
		while (i < STRINGLEN) {
			scanf("%d",&s);
			if ((string[i++] = s) == 0) break;
		}
		string[i] = '\0';
	}
}
drawstrings(x,y,string)
int x,y;
char string[];
{
	double magx(),magy();
	XDrawString(display,win,gc,(int)magx(x),(int)magy(y),string,(int)strlen(string));
}

set_misc(sv,sc,xfactor,yfactor,la,lav,arrowtail)
	int sv[],lav[],arrowtail[];
	double sc[],*xfactor,*yfactor,la[];
{
	*xfactor = ((double)(sv[1]-sv[0]))/(sc[1]-sc[0]);
	*yfactor = ((double)(sv[3]-sv[2]))/(sc[3]-sc[2]);
	la[0] = (sc[1]-sc[0])/20.0;
	la[1] = (sc[3]-sc[2])/20.0;
	lav[0] = 5;
	lav[1] = 5;
	arrowpoint[0][1].x =   (sv[1]-sv[0])/8.0;
	arrowpoint[1][1].y = - (sv[3]-sv[2])/8.0;
	arrowtail[0] = sv[0]+(sv[1]-sv[0])/2-(arrowpoint[0][1].x)/2;
	arrowtail[1] = sv[3]+80;
	arrowtail[2] = sv[0]-80;
	arrowtail[3] = sv[2]+(sv[3]-sv[2])/2-(arrowpoint[1][1].y)/2;
}

struct2data(xddp,sv,sc,la,lav,fi,kse,fontname, arrowtail,flag)
	struct xddp_dat xddp;
	int sv[],lav[],arrowtail[],*kse,*flag,fi[];
	double sc[],la[];
	char *fontname;
{
	int i;
	for (i = 0; i < 4; i++){ sv[i] = xddp.sv[i];}
	for (i = 0; i < 4; i++){ sc[i] = xddp.sc[i];}
	for (i = 0; i < 2; i++){ la[i] = xddp.la[i];}
	for (i = 0; i < 2; i++){ lav[i] = xddp.lav[i];}
	for (i = 0; i < 4; i++){ arrowtail[i] = xddp.arrowtail[i];}
	for (i = 0; i < 2; i++){ fi[i] = xddp.fi[i];}
	*kse = xddp.kse;
	*flag = xddp.flag;
	strcpy(fontname,xddp.font);
	font = XLoadQueryFont(display,fontname);
	XSetFont(display,gc,font->fid);
}

int no_wait_key()
{
	int ix,iy;
	int key;
	char string[10];
	Bool test;
	XEvent event;
	test = XCheckMaskEvent(display,KeyPressMask,&event);
	if (test == True &&
		event.type == KeyPress &&
			XLookupString(&event, string, 10, &key, NULL) == 1){
			xddpret.key = (int)string[0];
	} else  xddpret.key = -1;
}

int key_in()
{
	int ix,iy;
	int key;
	char string[10];
	Bool test;
	XKeyEvent event;
	XSelectInput(display,win,KeyPressMask);
	while (1) {
		XNextEvent(display,&event);
		if (event.type == KeyPress &&
			XLookupString(&event, string, 10, &key, NULL) == 1){
			xddpret.key = (int)string[0];
			break;
		}
	}
}


int bitmapout()
{
	Window root;
	Pixmap bmap,pc;
	int depth;
	int x, y;
	unsigned int width, height;
	unsigned b,d;
	GC gc0;

	depth = DefaultDepth(display, 0);
	printf("%d\n", depth);

	XGetGeometry(display, win, &root, &x, &y, &width, &height, &b, &d); 
    bmap = XCreatePixmap(display, win, width, height, 1);
	gc0 = XCreateGC(display, bmap, 0,0);
	switch (depth){
		case 1:
		    XCopyArea(display, win, bmap, gc0, 0,0,
				width, height, 0, 0);
    		XWriteBitmapFile(display, "xddp.xbm", bmap,
				width, height, 0,0);
			break;
		case 8:
		    XCopyPlane(display, win, bmap, gc0, 0,0,
				width, height, 0, 0, 1);
    		XWriteBitmapFile(display, "xddp.xbm", bmap,
				width, height, 0,0);
			break;
		default:
			printf("unknown depth %d\n", depth);
			break;
	}
    XFreePixmap(display, bmap);
    XFreeGC(display, gc0);
}

mkarrow()
{
	char string[10];
	int key;
	int i, button;
	int x0,y0,x1,y1;
	int ix0,iy0;
	int ix1,iy1;
	int ix,iy;
	double xa,ya,xaa, yaa, xb,yb;
	XPoint xp[4];
	double arg;
	double r = 25.0;
	double delta = 15.0 * M_PI / 180.0;
	double x1d,y1d,xaad, yaad,xad,yad,xbd,ybd;
	double ca,sa;
	XEvent event;
	XKeyEvent kevent;

	XSelectInput(display,win,ButtonPressMask|KeyPressMask);
	i = 0;
	while (1){
		XNextEvent(display,&event);
		if (event.type == ButtonPress) {
			i++;
			switch (i){
				case 1:
					ix0 = event.xbutton.x;
					iy0 = event.xbutton.y;
					printf("%d,%d\n",ix0,iy0);
					break;
				case 2:
					ix1 = event.xbutton.x;
					iy1 = event.xbutton.y;
					printf("%d,%d\n",ix1,iy1);
					break;
			}
		} else if (event.type == KeyPress &&
			XLookupString(&kevent, string, 10, &key, NULL) == 1){
			xddpret.key = (int)string[0];
			return(-3);
		} 
		if (i == 2) break;
	}
	ix = ix1 - ix0; iy = -(iy1 - iy0);
	arg = atan2((double)(ix),(double)(-iy));
	arg -= M_PI / 2.0; 
	sa = sin(arg); ca = cos(arg);
	x1d = (double)(ix) * ca + (double)(iy) * sa; 
	y1d = -(double)(ix) * sa + (double)(iy) * ca; 
	xad = x1d - r * cos(delta);
	yad = y1d + r * sin(delta);
	xaad = x1d - r * 0.5;
	yaad = y1d;
	xbd = xad;
	ybd = y1d - r * sin(delta);
	xa = xad * ca - yad * sa;
	ya = xad * sa + yad * ca;
	xaa = xaad * ca - yaad * sa;
	yaa = xaad * sa + yaad * ca;
	xb = xbd * ca - ybd * sa;
	yb = xbd * sa + ybd * ca;
	xp[0].x = (short int)(x1d * ca - y1d * sa)+ix0;
	xp[0].y = -(short int)(x1d * sa + y1d * ca)+iy0;
	xp[1].x = (short int)xa + ix0;
	xp[1].y = -(short int)ya + iy0;
	xp[2].x = (short int)xaa + ix0;
	xp[2].y = -(short int)yaa + iy0;
	xp[3].x = (short int)xb + ix0;
	xp[3].y = -(short int)yb + iy0;
	XFillPolygon(display, win, gc, xp, 4, Nonconvex, CoordModeOrigin);
	xddpret.key = -1;
	}
