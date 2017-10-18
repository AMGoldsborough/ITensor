include ../this_dir.mk
include ../options.mk
################################################################

TENSOR_HEADERS=core.h

#################################################################

#Mappings --------------
REL_TENSOR_HEADERS=$(patsubst %,$(ITENSOR_INCLUDEDIR)/%, $(TENSOR_HEADERS))

#Define Flags ----------
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: disordered_spinhalf disordered_spinhalfQ DMRG_dMERA DMRG_dMERA_IQ DMRG_dMERA_noee DMRG_dMERA_IQ_noee sdrg spinhalfssd_test

debug: dmrg-g  

all: disordered_spinhalf disordered_spinhalfQ DMRG_dMERA DMRG_dMERA_IQ DMRG_dMERA_noee DMRG_dMERA_IQ_noee sdrg spinhalfssd_test

disordered_spinhalf: disordered_spinhalf.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) disordered_spinhalf.o -o disordered_spinhalf $(LIBFLAGS)

disordered_spinhalfQ: disordered_spinhalfQ.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) disordered_spinhalfQ.o -o disordered_spinhalfQ $(LIBFLAGS)

DMRG_dMERA: DMRG_dMERA.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) DMRG_dMERA.o -o DMRG_dMERA $(LIBFLAGS)

DMRG_dMERA_IQ: DMRG_dMERA_IQ.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) DMRG_dMERA_IQ.o -o DMRG_dMERA_IQ $(LIBFLAGS)

DMRG_dMERA_noee: DMRG_dMERA_noee.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) DMRG_dMERA_noee.o -o DMRG_dMERA_noee $(LIBFLAGS)

DMRG_dMERA_IQ_noee: DMRG_dMERA_IQ_noee.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) DMRG_dMERA_IQ_noee.o -o DMRG_dMERA_IQ_noee $(LIBFLAGS)

sdrg: sdrg.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) sdrg.o -o sdrg $(LIBFLAGS)

spinhalfssd_test: spinhalfssd_testsdrg.o $(ITENSOR_LIBS) $(REL_TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) spinhalfssd_test.o -o spinhalfssd_test $(LIBFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	disordered_spinhalf disordered_spinhalfQ DMRG_dMERA DMRG_dMERA_IQ DMRG_dMERA_noee DMRG_dMERA_IQ_noee sdrg spinhalfssd_test
