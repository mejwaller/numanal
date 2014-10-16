#compiler
CC=g++

#includes
INCPATH += -I. -I/usr/local/qt/include 

#compilation options
CFLAGS+=-g -Wall -ansi
CPPFLAGS+=${CFLAGS}

#link options
#LDFLAGS+=

#link libraries
LIBS+= /usr/local/qt/lib/libqt.so ./numlib.a ./qtutils.a

.SUFFIXES:	.cpp

.cpp.o:
	$(CC) $(CPPFLAGS) $(INCPATH) -c $<


numanal: main.o numlib.a qtutils.a
	$(CC) -o numanal main.o $(LIBS)
 
numlib.a: numlib.a(Bairstow1.o) \
 numlib.a(Bisection1.o) numlib.a(FalsePos1.o) \
 numlib.a(Horner1.o) numlib.a(Matrix1.o) numlib.a(Muller1.o) \
 numlib.a(Newton1.o) numlib.a(Secant1.o) numlib.a(Interpolate.o) \
 numlib.a(MathUtils.o) numlib.a(Bezier1.o)
 
qtutils.a: qtutils.a(ChartWidget1.o)

main.o: main.cpp Bisection1.hpp BaseUnaryNonLinearSolver.hpp Secant1.hpp \
  FalsePos1.hpp Newton1.hpp Muller1.hpp Horner1.hpp Bairstow1.hpp \
  Matrix1.hpp Interpolate.hpp ChartWidget1.hpp Bezier1.hpp MathUtils.hpp

Bairstow1.o: Bairstow1.cpp Bairstow1.hpp BaseUnaryNonLinearSolver.hpp
 
Bisection1.o: Bisection1.cpp Bisection1.hpp BaseUnaryNonLinearSolver.hpp

FalsePos1.o: FalsePos1.cpp FalsePos1.hpp BaseUnaryNonLinearSolver.hpp

Horner1.o: Horner1.cpp Horner1.hpp

Matrix1.o: Matrix1.cpp Matrix1.hpp BaseUnaryNonLinearSolver.hpp

Muller1.o: Muller1.cpp Muller1.hpp BaseUnaryNonLinearSolver.hpp

Newton1.o: Newton1.cpp Newton1.hpp BaseUnaryNonLinearSolver.hpp

Secant1.o: Secant1.cpp Secant1.hpp BaseUnaryNonLinearSolver.hpp

Interpolate.o: Interpolate.cpp Interpolate.hpp Matrix1.hpp

MathUtils.o: MathUtils.cpp MathUtils.hpp

Bezier1.o: Bezier1.cpp Bezier1.hpp MathUtils.hpp

ChartWidget1.o: ChartWidget1.cpp ChartWidget1.hpp

