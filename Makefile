###################################################################
# This Makefile was created using the bat-project script
# for project fatjetModel
# bat-project is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CXXFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# List of all class (models) sources used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS = fatjetModel.cxx

# List of all program sources used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
PRGSRCS = runfatjetModel.cxx get_prediction.cxx #make_limits.cxx

# compiler and flags
CXX       = g++
CXXFLAGS  = -g -O2 -Wall -fPIC -Wno-deprecated 
LD        = /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/x86_64-pc-linux-gnu/bin/ld
LDFLAGS   = -g -O2 

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

CXXFLAGS += $(shell bat-config --cflags)
LIBS := $(shell bat-config --libs)

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS = $(addsuffix .o,$(basename $(CXXSRCS)))
MYPROGS = $(basename $(PRGSRCS))
PRGOBJS = $(addsuffix .o,$(basename $(PRGSRCS)))

GARBAGE = $(CXXOBJS) $(PRGOBJS) link.d $(MYPROGS)

# targets
all : $(MYPROGS)

.PHONY : all clean print

link.d : $(addsuffix .h,$(basename $(CXXSRCS))) $(CXXSRCS) $(PRGSRCS)
	$(CXX) -MM $(CXXFLAGS) $(filter-out %.h,$^) > link.d;
	@$(foreach prog,$(MYPROGS), echo $(prog) : $(prog).o >> link.d;)

-include link.d

$(CXXOBJS) $(PRGOBJS) :
	$(CXX) $(CXXFLAGS) -c $(filter $(basename $@).%,$(filter-out %.h,$^)) -o $@

$(MYPROGS) : $(CXXOBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(GARBAGE)

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS) $(PRGSRCS)
	@echo c++ objs  : $(CXXOBJS) $(PRGOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)
