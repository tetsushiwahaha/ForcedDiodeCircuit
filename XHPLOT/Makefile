# Makefile: xhplot for Solaris
#
CC = clang
CFLAGS = -O -I$(INCLUDE_DIR)
#AR = /usr/xpg4/bin/ar
#AR = /usr/bin/ar
#AR = /usr/ccs/bin/ar
AR = /opt/local/bin/ar
INCLUDE_DIR = /usr/X11R6/include
LIBNAME = libxddp
OBJ = xhplot

install: $(OBJ).o
	$(AR) r $(LIBNAME).a $(OBJ).o
	mv $(LIBNAME).a $(HOME)/lib/$(LIBNAME).a

clean: 
	\rm *.o *.a
