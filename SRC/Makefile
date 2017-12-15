


###########################################################################################
# COMPILERS OPTIONS
CXX=g++ 
CXX:= ccache $(CXX)
CXXFLAGS = -O3 -std=c++11 -g #-v 

###########################################################################################
BINDIR:=bin
SRCDIR:=src
BUILDDIR:=build

###########################################################################################
# FILES 
EXT=cc
SOURCES:=$(wildcard $(SRCDIR)/*.$(EXT))
HEADERS:=$(wildcard $(SRCDIR)/*.h)
OBJECTS:=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(EXT)=.o))


###########################################################################################
# -L: FOLDER LIBRARY
#
LFLAGS+=-L./lib  
LFLAGS+=-L/usr/lib

# -l: LIBRARY 
lLIBS+= -llog4cxx
lLIBS+= -larmadillo


###########################################################################################
# -I: INCLUDES 
IFLAGS+=-I./include
IFLAGS+=-I./src


###########################################################################################
# CMD OPTIONS 

#sufix_1="NOPARALLEL"
#sufix_2="NOOPTIMIZE"
#sufix_3="NODEBUG"
#
#ifneq ("$(PARALLEL)","NOPARALLEL")
#CXXFLAGS+=-D PARALLEL 
#sufix_1="PARALLEL"
#endif
#
#ifneq ("$(OPTIMIZE)","NOOPTIMIZE") 
#CXXFLAGS+=-D OPTIMIZE
#sufix_2="OPTIMIZE"
#endif
#
#ifeq ("$(DEBUG)","DEBUG") 
#CXXFLAGS+= -g -Wall
#CXXFLAGS+=-D DEBUG
#sufix_3="DEBUG"
#endif
#



###########################################################################################
###########################################################################################
TARGET=all ################################################################################

# make
$(TARGET): $(OBJECTS)
	$(info LINKING: $@) 
	@$(CXX) -o $(BINDIR)/OWR $(OBJECTS) $(LFLAGS) $(lLIBS) $(SFLAGS) $(HLDFLAGS)  
	#@$(CXX) -o $(BINDIR)/OWR$(sufix_1)_$(sufix_2)_$(sufix_3) $(OBJECTS) $(LFLAGS) $(lLIBS) $(SFLAGS) $(HLDFLAGS)  

###########################################################################################
$(BUILDDIR)/%.o : $(SRCDIR)/%.$(EXT)
	$(info compiling: $@) 
	@$(CXX) $(CXXFLAGS) -c $<  -o $@ $(HCPPFLAGS) $(IFLAGS)  


###########################################################################################
.PHONY: clean
clean:
	@$(RM) -f $(TARGET) $(OBJECTS) 2>/dev/null || true
	@rm $(BINDIR)/aw_$(sufix_1)_$(sufix_2)_$(sufix_3)  2>/dev/null || true
	clear

###########################################################################################
list:
	@echo
	$(info LIST OF SOURCES  + HEADERS + OBJECTS:) 
	@echo
	@echo $(SOURCES)
	@echo
	@echo $(HEADERS)
	@echo
	@echo $(OBJECTS)
	@echo
###########################################################################################
docu:
	@echo
	$(info Generating doxygen docu in : /doc/DOXYGEN_DOCUMENTATION) 
	@echo
	cd src; doxygen ../doc/config-file-doxygen; mv DOXYGEN_DOCUMENTATION ../doc

