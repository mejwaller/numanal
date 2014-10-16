#compiler
CC=gcc

#includes
INCLUDE+=.

#compilation options
CFLAGS+=-g -Wall -ansi
CPPFLAGS+=${CFLAGS}

#link options
#LDFLAGS+=

#link libraries
#LDLIBS+=

nr: cmain.o nr.o matrices.o
	$(CC) -o nr cmain.o nr.o matrices.o

cmain.o: cmain.c matrices.h nr.h

nr.o: nr.c nr.h

matrices.o: matrices.c matrices.h nr.h


