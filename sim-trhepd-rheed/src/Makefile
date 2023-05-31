RM = rm
ifeq ($(OS),Windows_NT)
    RM = cmd.exe /C del
endif

BLAS_CFLAGS=
BLAS_LIBS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

FC = gfortran
FCFLAGS=-O3 -g -Wall -Warray-temporaries -fmax-stack-var-size=1000 -fopenmp $(BLAS_CFLAGS)
LDFLAGS=$(BLAS_LIBS) -lgfortran
#FCFLAGS=-Wall -pedantic -fbounds-check -g -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
BULK = bulkm.o bulkio.o bulk.o bulksub.o
SURF = surfm.o surflf.o surfio.o surfsub.o domainsum.o surf_prkn.o matcomp.o surf_magnus.o matgs.o surf_rk.o h211b_control.o 
DOMAIN = domain.o domainsum.o
XYZ = xyzb.o strfac.o asf.o 
SOURCES = asf.o strfac.o scpot.o trmatg.o gcmi.o 
POTCALC = U0.o bulkU0.o surfU0.o bulksub.o surfsub.o asf.o strfac.o scpotU0.o
LAPACK = 
# LAPACK = -llapack -lblas

ALL: bulk surf xyz potcalc

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<
%.o: %.f
	$(FC) $(FCFLAGS) -o $@ -c $<

surf_prkn.o: matcomp.o h211b_control.o surfio.o
surf_magnus.o: matcomp.o h211b_control.o surfio.o
surf_rk.o : matcomp.o surfio.o matgs.o
surflf.o : surf_prkn.o surf_magnus.o matcomp.o surfio.o surf_rk.o

bulk: $(BULK) $(SOURCES) $(LAPACK)
	$(FC) $(FCFLAGS) -o bulk.exe $(BULK) $(SOURCES) $(LAPACK) $(LDFLAGS)

surf: $(SURF) $(SOURCES) $(LAPACK)
	$(FC) $(FCFLAGS) -o surf.exe $(SURF) $(SOURCES) $(LAPACK) $(LDFLAGS)

domain: $(DOMAIN) $(LAPACK)
	$(FC) $(FCFLAGS) -o domain.exe $(DOMAIN) $(LDFLAGS)

xyz: $(XYZ) $(LAPACK)
	$(FC) $(FCFLAGS) -o xyz.exe $(XYZ) $(LDFLAGS)

potcalc: $(POTCALC) $(LAPACK)
	$(FC) $(FCFLAGS) -o potcalc.exe $(POTCALC) $(LDFLAGS)

clean:
	-$(RM) *.exe *.o *.mod

