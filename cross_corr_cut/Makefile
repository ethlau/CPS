OBJ1 = halo_model.o mf.o\
	nrutil.o qromb.o polint.o trapzd.o spline.o\
        auxiliary_sin_cos_integrals_fi_gi.o\
        power_series_Cin.o power_series_Si.o\
        sin_cos_integrals.o
EXE1 = halo_model

#linker
#LD = g++ 
LD = icpc

#c compiler
CC =icc

#cpp compiler
CPP =icpc

#CPP compiler options
CPPOPT = -O3 #

#C compiler options
COPT = -O3  -s # -bmaxstack:512000000 -g

#f90 compiler
#FC90 = xlf

#linking
GSL = -lm -lgsl -lgslcblas

ALL = $(EXE1) $(EXE2)

#compilation
.C.o:
	$(CC) $(COPT) $(FCFLAGS) $(CPPFLAGS) -c $.C

%.o: %.C
	$(CPP) $(CPPOPT) -c $*.C

all: $(ALL)

$(EXE1): $(OBJ1)
	$(LD) $(CPPOPT) -o $(EXE1) $(OBJ1) $(GSL)
	echo "   Done"

.FAILED:
.DONE:

clean:
	rm -f $(ALL) *.o *~

#dependencies

