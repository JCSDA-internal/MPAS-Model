.SUFFIXES: .F .o

CHEMISTRY=-DDO_CHEMISTRY
#CHEMISTRY=
PHYSICS=-DDO_PHYSICS
#PHYSICS=

OBJS = mpas_atm_core.o \
       mpas_atm_core_interface.o \
       mpas_atm_dimensions.o \
       mpas_atm_threading.o

all: physcore_init chemcore physcore dycore diagcore atmcore utilities

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	( cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.atmosphere in_defaults=true )
	( cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.atmosphere stream_list.atmosphere. listed )

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) < ../Registry_processed.xml )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

physcore_init:
	(cd physics;  $(MAKE) core_physics_init )

physcore: mpas_atm_dimensions.o
	( cd physics; $(MAKE) all )
	( mkdir libphys; cd libphys; ar -x ../physics/libphys.a )
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*TBL .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*DATA* .)

chemcore: physcore_init
	( cd chemistry; $(MAKE) all )
	( mkdir libchem; cd libchem; ar -x ../chemistry/libchem.a )

dycore: mpas_atm_dimensions.o chemcore physcore
	( cd dynamics; $(MAKE) all PHYSICS="$(PHYSICS)" CHEMISTRY="$(CHEMISTRY)" )

diagcore: chemcore physcore dycore
	( cd diagnostics; $(MAKE) all )

utilities: physcore
	( cd utils; $(MAKE) all )

atmcore: chemcore physcore dycore diagcore $(OBJS)
	ar -ru libdycore.a $(OBJS) dynamics/*.o libphys/*.o libchem/*.o diagnostics/*.o

mpas_atm_core_interface.o: mpas_atm_core.o

mpas_atm_core.o: dycore diagcore mpas_atm_threading.o

mpas_atm_dimensions.o:

clean:
	( cd physics; $(MAKE) clean )
	( cd chemistry; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	( cd diagnostics; $(MAKE) clean )
	( cd utils; $(MAKE) clean )
	( cd ../..; rm -f *TBL )
	( cd ../..; rm -f *DATA* )
	$(RM) -r libphys
	$(RM) -r libchem
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CHEMISTRY) $(PHYSICS) $(CPPINCLUDES) -I./inc $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./chemistry -I./physics -I./dynamics -I./diagnostics -I./chemistry/chemistry_opp -I./physics/physics_wrf -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(CHEMISTRY) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I./inc -I../framework -I../operators -I./chemistry -I./physics -I./dynamics -I./diagnostics -I./chemistry/chemistry_opp -I./physics/physics_wrf -I../external/esmf_time_f90
endif
