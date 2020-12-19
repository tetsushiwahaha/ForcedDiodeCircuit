/*
	commands 
*/
#define	INIT_DRAW 0
#define	REDRAW 1
#define	PSET 2
#define	PRELINE 3
#define	LINE 4
#define	CLS 5
#define	LINEATT 6
#define	ARCSTRUCT 7
#define	ARC 8
#define	INIT_SETP 10
#define	SETPOINTS 11
#define	PSETS 12
#define	FILL 13
#define	LINES 14
#define	BITMAP 15
#define	FILLOpaque 16
#define	CLIPMASK 17
#define	CLIPFULL 18
#define	ARROW 19
#define	SETARROW 20
#define	BITMAPS 21
#define	SETBITMAP 22
#define	WINDOW 23
#define	INIT_DRAW_RV 25
#define SETFONT 26
#define STRING 27
#define SETSTRING 28
#define STRINGS 29
#define FRAME 30
#define DIRECT 31
#define CYLINDER 32
#define CYLPSET 33
#define CYLLINE 34
#define CYLFRAME 35
#define DIRECT_RV 36
#define	CLOSE 99
#define	WAIT 999
#define KEYIN 9996
#define STATICKEYIN 9995
#define	POINTER 9999
#define	CYLPOINT 9997
#define	POINT 9998
#define	MAGNIFY 1203
#define BITMAPOUT 500
#define MKARROW 501
/*
	set xhplot.init
*/
#define	AXES 1
#define	GRID 2
#define	ABNORMAL 3
#define	AXES_U 4
#define	AXES_R 5
#define	AXES_UR 6
#define	AXES_C 7
#define	FRAMEONLY 8
#define	OKDRAW 1
#define	NOTDRAW 2
#define	SETAGAIN 3
#define	NOLINE 4
#define	NOMOJI 5
/*
	line style
*/
#define	SOLID 1
#define	DASH 2
#define	DOTTED 3
#define	DOT_DASHED 4
#define	SHORT_DASHED 5
#define	LONG_DASHED 6
#define	ODD_DASHED 7
/*
	colors avairable
*/
#define	BLACK 0
#define	BLUE 1
#define	RED 2
#define	MAGENTA 3
#define	GREEN 4
#define	CYAN 5
#define	YELLOW 6
#define	WHITE 7
#define	STEELBLUE 8
#define	NAVYBLUE 9
#define	NAVY 10
#define	INDIANRED 11
#define	PINK 12
#define	DARKGREEN 13
#define	SPRINGGREEN 14
#define	TURQUOISE 15
/*
	arrow
*/
#define	HORIZONTAL 0
#define	VERTICAL 1
#define	TAIL 0
#define	HEAD 1
#define	NORMAL 0
#define	INVERSE 1

struct xddp_dat {
	int line_wid,line_att;
	int sv[4];
	double sc[4];
    double la[2];
    int lav[2];
    int fi[2];
    int kse;
    char font[100];
    int arrowtail[4];
    int flag;
} ;

struct xddp_arc{
	int height;
	int width;
	int angle1;
	int angle2;
	int icol;
};


struct xddp_return {
	double x,y;
	int key;
};

struct xddp_dat xddp;
struct xddp_return xddpret;
struct xddp_arc xddparc;

