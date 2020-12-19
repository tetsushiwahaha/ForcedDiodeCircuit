#include<sys/time.h>
#include<stdio.h>
#include<sys/termio.h>
#define STDIN 0
suzu(ret,c)
int *ret;
char *c;
{
	int readfd;
	struct timeval timeout;
	struct termio stdio;
	struct termio status;

	ioctl(STDIN,TCGETA,&stdio);
	ioctl(STDIN,TCGETA,&status);
	stdio.c_lflag &= ~ICANON;
	stdio.c_cc[VEOF] = 1;
	ioctl(STDIN,TCSETA,&stdio);

	timeout.tv_sec = 0;
	timeout.tv_usec = 0;
	readfd = 1;
	*ret = select(1,&readfd,0,0,&timeout);
	if(*ret>0) read(STDIN,c,1);
	ioctl(STDIN,TCSETA,&status);
}
