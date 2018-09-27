#linker
#LD = g++ 
LD = icpc

#CC = gcc
#CXX= g++
#FC = gfortran

CC = icc
CXX= icpc
FC = ifort

#
#CPP compiler options
CPPOPT = -O3 -I./DK15 -std=c++11

#C compiler options
COPT = -O3  -s # -bmaxstack:512000000 -g

#f90 compiler
#FC90 = xlf

#linking
GSL = -lm -lgsl -lgslcblas
#compilation

DK15_SRCS = halo_conc_DK15.cpp nrutil.cpp polint.cpp qromb.cpp spline.cpp trapzd.cpp
DK15_OBJS = $(patsubst %.cpp,DK15/%.o,$(DK15_SRCS))

DK15/%.o: DK15/%.cpp DK15/%.h
	$(CXX) -O2 -I. -I./DK15 -c $< -o $@


objects= fftlog.o cdgamma.o drfftb.o drfftf.o drffti.o
save_yl: save_yl.cpp $(objects) gas_model.o $(DK15_OBJS)
	$(CXX) $(CPPOPT) -o save_yl save_yl.cpp $(objects) gas_model.o $(DK15_OBJS) $(GSL) -lifcore -lifport

gas_model.o: gas_model.cpp gas_model.h
	$(CXX) $(CXXFLAGS) -c gas_model.cpp

	echo "   Done"

.FAILED:
.DONE:

clean:
	rm -f *.o DK15/*.o save_yl *~

#dependencies

