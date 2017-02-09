CXX	= g++
FORTRAN	= gfortran -fno-automatic
RC1	= `root-config --cflags`
RC2	= `root-config --cflags --glibs`
INC	= $(CLHEP_BASE_DIR)/include
LIBPATH	= $(CLHEP_BASE_DIR)/lib
LIB	= -lCLHEP -lGeom -lGeomPainter -lGed -lgfortran

all: monopoles

clean:
	rm -f *.o

d2_integral.o: d2_integral.f
	$(FORTRAN) -o d2_integral.o -c d2_integral.f
gaisser.o: gaisser.f
	$(FORTRAN) -o gaisser.o -c gaisser.f
hrndg2.o: hrndg2.f
	$(FORTRAN) -o hrndg2.o -c hrndg2.f
monopoles.o: monopoles.cxx
	$(CXX) -o monopoles.o -c $(RC1) -I$(INC) monopoles.cxx
detector.o: detector.cxx
	$(CXX) -o detector.o -c $(RC1) -I$(INC) detector.cxx
monopoles: monopoles.o detector.o d2_integral.o gaisser.o hrndg2.o 
	$(CXX) monopoles.o detector.o d2_integral.o gaisser.o hrndg2.o $(RC2) -L$(LIBPATH) $(LIB) -o monopoles -O3

