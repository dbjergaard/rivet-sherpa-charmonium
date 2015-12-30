CC=g++

LIBDIR:=$(shell rivet-config --libdir)

INCDIR=$(PWD)/include
RIVETINCDIR:=$(shell rivet-config --cppflags)
LDFLAGS:=$(shell rivet-config --ldflags) -ldl
WFLAGS= -Wall -Wno-long-long -Wno-format #-Werror=uninitialized -Werror=delete-non-virtual-dtor  -Wno-unused-local-typedefs
CFLAGS= -g3 -fno-inline -I$(INCDIR) $(RIVETINCDIR) -pedantic -ansi $(WFLAGS) -O0 -Wl,--no-as-needed -lRivet

.PHONY: all
# Code related rules
all: rivet-lib #Pythia8/pythia8
rivet-lib: RivetMC_GENSTUDY_CHARMONIUM.so RivetMC_GENSTUDY_JPSI_VEC.so libBOOSTFastJets.so
Pythia8/pythia8: Pythia8/pythia8.o
	$(CC) $^ -o $@ $(LDFLAGS) -lpythia8 -lHepMC 
Pythia8/pythia8185: Pythia8/pythia8185.o
	$(CC) $^ -o $@ $(LDFLAGS) -lpythia8 -lpythia8tohepmc -lHepMC -L ${HOME}/rivet/local/Pythia8185/lib -L ${HOME}/rivet/local/Pythia8185/lib/archive -llhapdfdummy
# Pythia8/pythia8210: Pythia8/pythia8210.o
# 	$(CC) $^ -o $@ $(LDFLAGS) -lpythia8 -lHepMC -L ${HOME}/rivet/local/Pythia8210/lib
# Pythia8/pythia8185.o: Pythia8/runPythia.C
# 	$(CC) $(CFLAGS) -I${HOME}/rivet/local/Pythia8185/include -c $^ -o $@
# Pythia8/pythia8210.o: Pythia8/runPythia.C
# 	$(CC) $(CFLAGS) -D__PYTHIA_8210__ -I${HOME}/rivet/local/Pythia8210/include -c $^ -o $@
Pythia8/pythia8.o: Pythia8/runPythia.C
	$(CC) $(CFLAGS) -D__PYTHIA_8210__ -c $^ -o $@
RivetMC_GENSTUDY_CHARMONIUM.so:  MC_GENSTUDY_CHARMONIUM.cc libBOOSTFastJets.so
	$(CC) -o "$@" -shared -fPIC $(CFLAGS) $< -lBOOSTFastJets $(LDFLAGS) -lfastjetcontribfragile -L ./
RivetMC_GENSTUDY_JPSI_VEC.so:  MC_GENSTUDY_JPSI_VEC.cc libBOOSTFastJets.so
	$(CC) -o "$@" -shared -fPIC $(CFLAGS) $< -lBOOSTFastJets $(LDFLAGS) -lfastjetcontribfragile -L ./
libBOOSTFastJets.so: src/BOOSTFastJets.cxx
	$(CC)  -shared -fPIC $(CFLAGS) $< -o $@ -lfastjet -lfastjettools $(LDFLAGS)
# Cleaning and installing rules
install:
	cp *.so $(LIBDIR)
clean:
	-rm -f *.o  *.so 
