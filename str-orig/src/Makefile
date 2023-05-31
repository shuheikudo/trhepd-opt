RM = rm
ifeq ($(OS),Windows_NT)
    RM = cmd.exe /C del
endif

BLAS_CFLAGS=
BLAS_LIBS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

FC = gfortran
FCFLAGS=-O3 -g -Wall -Warray-temporaries -fmax-stack-var-size=1000 $(BLAS_CFLAGS)
LDFLAGS=$(BLAS_LIBS) -lgfortran
#FCFLAGS=-Wall -pedantic -fbounds-check -g -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
BULK = bulkm.f90 bulkio.f90 bulk.f90 bulksub.f90
SURF = surfm.f90 surf.f90 surfio.f90 surfsub.f90 domainsum.f90
DOMAIN = domain.f90 domainsum.f90
XYZ = xyzb.f90 strfac.f90 asf.f90 
SOURCES = asf.f90 strfac.f90 scpot.f90 trmatg.f90 gcmi.f90
POTCALC = U0.f90 bulkU0.f90 surfU0.f90 bulksub.f90 surfsub.f90 asf.f90 strfac.f90 scpotU0.f90
LAPACK = 
# LAPACK = -llapack -lblas

ALL: bulk surf xyz potcalc

bulk: $(BULK) $(SOURCES)
	$(FC) $(FCFLAGS) -o bulk.exe $(BULK) $(SOURCES) $(LAPACK) $(LDFLAGS)

surf: $(SURF) $(SOURCES)
	$(FC) $(FCFLAGS) -o surf.exe $(SURF) $(SOURCES) $(LAPACK) $(LDFLAGS)

domain: $(DOMAIN)
	$(FC) $(FCFLAGS) -o domain.exe $(DOMAIN) $(LDFLAGS)

xyz: $(XYZ)
	$(FC) $(FCFLAGS) -o xyz.exe $(XYZ) $(LDFLAGS)

potcalc: $(POTCALC)
	$(FC) $(FCFLAGS) -o potcalc.exe $(POTCALC) $(LDFLAGS)

clean:
	-$(RM) *.exe 

