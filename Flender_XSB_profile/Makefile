#CC = /opt/local/bin/gcc-mp-4.8  # gcc
#CXX= /opt/local/bin/g++-mp-4.8  # g++
#FC = gfortran
CC = icc
CXX= icpc
FC = ifort
#CC=gcc
#CXX=g++
#FC=gfortran

CFLAGS = -O2 
CXXFLAGS = -O2 -std=c++11

#CXXFLAGS = -O0 -g
#CXXFLAGS += -DMPI_PARALLEL
#CXXFLAGS += -DOUTPUT_NFW_DENSITY
#CXXFLAGS += -DLONGIDS
#CXXFLAGS += -DEXCLUDE_SUBHALOS
#CXXFLAGS += -DMASS_SELECTION=1e14
#CXXFLAGS += -DROCKSTAR_CONCENTRATION
#FLINE = -D'LINE_FITS_FILE="/Volumes/razoku01/xc_RT2048/code/Xraycalc/atomdb_v3.0.9/apec_line.fits"'
#FCOCO = -D'COCO_FITS_FILE="/Volumes/razoku01/xc_RT2048/code/Xraycalc/atomdb_v3.0.9/apec_coco.fits"'
FLINE = -D'LINE_FITS_FILE="/home/fas/nagai/etl28/programs/Xrays/atomdb/atomdb_v3.0.9/apec_line.fits"'
FCOCO = -D'COCO_FITS_FILE="/home/fas/nagai/etl28/programs/Xrays/atomdb/atomdb_v3.0.9/apec_coco.fits"'
#FLINE = -D'LINE_FITS_FILE="/Users/ethlau/Research/Xrays/atomdb_v3.0.9/apec_line.fits"'
#FCOCO = -D'COCO_FITS_FILE="/Users/ethlau/Research/Xrays/atomdb_v3.0.9/apec_coco.fits"'

#CFLAGS += -I/Users/masato/gsl-2.4/include -I/usr/local/include
#CXXFLAGS += -I/Users/masato/gsl-2.4/include -I/usr/local/include -I./DK15
#CLIBS = -L/Users/masato/gsl-2.4/lib -lgsl -lgslcblas -lm -L/usr/local/lib -lcfitsio -lfftw3
#CLIBS += -lgfortran /opt/local/lib/gcc48/libquadmath.a

CFLAGS += -I/home/fas/nagai/etl28/local/include
CXXFLAGS += -I/home/fas/nagai/etl28/local/include -I./DK15
CLIBS = -lgsl -lgslcblas -lm 
CLIBS += -L/home/fas/nagai/etl28/local/lib/ -lfftw3 -lcfitsio -lifcore -lifport

#CFLAGS += -I/usr/local/include
#CXXFLAGS += -I/usr/local/include -I./DK15
#CLIBS = -lgsl -lgslcblas -lm
#CLIBS += -lfftw3 -lcfitsio -lgfortran

APEC_SRCS = Apec.c atomdb_make_spectrum.c calc_continuum.c calc_lines.c messages.c readapec.c read_continuum_data.c read_fits_spectrum.c read_line_data.c read_parameters.c gaussianLine.c
APEC_OBJS = $(patsubst %.c,Apec/%.o,$(APEC_SRCS))

DK15_SRCS = halo_conc_DK15.cpp nrutil.cpp polint.cpp qromb.cpp spline.cpp trapzd.cpp
DK15_OBJS = $(patsubst %.cpp,DK15/%.o,$(DK15_SRCS))

DK15/%.o: DK15/%.cpp DK15/%.h
	$(CXX) -O2 -I. -I./DK15 -c $< -o $@

Apec/%.o: Apec/%.c Apec/%.h
	$(CC) $(CFLAGS) $(FCOCO) $(FLINE) -I. -I./Apec -c $< -o $@

objects= fftlog.o cdgamma.o drfftb.o drfftf.o drffti.o
save_xl: save_xl.cpp $(objects) gas_model.o xray.o $(APEC_OBJS) $(DK15_OBJS)
	$(CXX) $(CXXFLAGS) -o save_xl save_xl.cpp $(objects) gas_model.o xray.o $(APEC_OBJS) $(DK15_OBJS) $(CLIBS)

xray.o: xray.c xray.h
	$(CC) $(CFLAGS) -c xray.c

gas_model.o: gas_model.cpp gas_model.h
	$(CXX) $(CXXFLAGS) -c gas_model.cpp


clean:
	/bin/rm -f *.o DK15/*.o Apec/*.o save_xl *~
