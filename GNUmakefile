 #############################################################################
 # Project: RooFit                                                           #
 # Package: RooFitCore                                                       #
 #    File: $Id: GNUmakefile.standalone,v 1.6 2004/04/05 22:38:34 wverkerke Exp $
 # Authors:                                                                  #
 #   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       #
 #   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 #
 #                                                                           #
 # Copyright (c) 2000-2004, Regents of the University of California          #
 #                          and Stanford University. All rights reserved.    #
 #                                                                           #
 # Redistribution and use in source and binary forms,                        #
 # with or without modification, are permitted according to the terms        #
 # listed in LICENSE (http://roofit.sourceforge.net/license.txt)             #
 #                                                                           #
 #   Standalone Makefile for RooFitModels                                    #
 #   ----------------------------------                                      #
 #                                                                           #
 #   Instructions                                                            #
 #       - Build RooFitCore first                                            #
 #         The RooFitCore subdir must be under the same directory            #
 #         as the RooFitModels subdir                                        #
 #                                                                           #
 #       - Review 'external configuration' section below                     #
 #         to match systems compilers setup                                  #
 #                                                                           #
 #       - Make sure the ROOTSYS environment variable is set and points      #
 #         to your ROOT release (3.02-07 or higher)                          #
 #                                                                           #
 #       - run 'gmake -f GNUMakefile.standalone <target>'                    #
 #         from the RooFitModels/ directory                                  #
 #                                                                           #
 #   Build targets                                                           #
 #     lib   - make libRooFitModels.a                                        #
 #     shlib - make libRooFitModels.so                                       #
 #     clean - delete all intermediate and final build objects               #
 #                                                                           #
 #   <NOTE TO USERS: This makefile is still a work in progress>              #
 #                                                                           #
 #   Copyright (C) 2002 University of California                             #
 #                                                                           #
 ############################################################################# 


# --- External configuration ---------------------------------
CC         = g++
CCFLAGS    = -O2 -Wno-deprecated -fPIC
MFLAGS     = -MM -Wno-deprecated
INCLUDES   = 
WORKDIR    = tmp
LIBDIR     = $(WORKDIR)
# -------------------------------------------------------------


# Internal configuration
PACKAGE=RooJpsiJpsiFit
LD_LIBRARY_PATH:=$(ROOTSYS)/lib
OBJDIR=$(WORKDIR)/objects
INCDIR=$(WORKDIR)/$(PACKAGE)
VPATH=$(INCDIR) $(OBJDIR)

INCLUDES += -I.. -I$(WORKDIR)/ -I$(ROOTSYS)/include -I$(PWD)/../RooFitCore/$(WORKDIR)
ROOTSYS  ?= ERROR_RootSysIsNotDefined
RDLLIST   = $(filter-out $(PACKAGE)_LinkDef.rdl,$(wildcard *.rdl))
CINTFILE  = $(PACKAGE)Cint.cc
CINTOBJ   = $(PACKAGE)Cint.o
LIBFILE   = $(LIBDIR)/lib$(PACKAGE).a
SHLIBFILE = $(LIBDIR)/lib$(PACKAGE).so

default: shlib

# List of all includes files to copy from rdl
HHLIST=$(patsubst %.rdl,%.hh,$(RDLLIST))

# List of all object files to build
OLIST=$(patsubst %.cc,%.o,$(wildcard *.cc))

# List of all dependency file to make
DLIST=$(patsubst %.rdl,%.d,$(RDLLIST))

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(INCDIR)/%.d: %.cc $(HHLIST)
	@echo "Making $@"
	@set -e; $(CC) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $< \
	          | sed 's/\($(notdir $*)\)\.o[ :]*/\1.o $(notdir $@) : /g' > $@; \
	        [ -s $@ ] || rm -f $@

# Implicit rule copying all RDL to INCDIR/HH
%.hh: %.rdl 
	@mkdir -p $(INCDIR)
	@cp $< $(INCDIR)/$@

# Implicit rule to compile all classes
%.o : %.cc 
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $(OBJDIR)/$(notdir $@) $(INCLUDES)


# Rule to make ROOTCINT output file
$(OBJDIR)/$(CINTOBJ): $(RDLLIST) 
	@mkdir -p $(INCDIR)
	@mkdir -p $(OBJDIR)
	@echo "Running rootcint"
	@cd $(INCDIR) ; $(ROOTSYS)/bin/rootcint -f $(CINTFILE) -c $(INCLUDES) $(HHLIST)
	@echo "Compiling $(CINTFILE)"
	@$(CC) $(CCFLAGS) -c $(INCDIR)/$(CINTFILE) -o $(OBJDIR)/$(CINTOBJ) $(INCLUDES)

# Rule to combine objects into a library
$(LIBFILE): $(OLIST) $(INCDIR)/$(CINTFILE)) $(patsubst %.cc,%.o,$(OBJDIR)/$(CINTFILE))
	@echo "Making $(LIBFILE)"
	@rm -f $(LIBFILE)
	@ar q $(LIBFILE) $(addprefix $(OBJDIR)/,$(OLIST) $(patsubst %.cc,%.o,$(CINTFILE))) 
	@ranlib $(LIBFILE)

# Rule to combine objects into a shared library
$(SHLIBFILE): $(OLIST) $(patsubst %.cc,%.o,$(OBJDIR)/$(CINTFILE))
	@echo "Making $(SHLIBFILE)"
	@rm -f $(SHLIBFILE)
	@$(CC) $(addprefix $(OBJDIR)/,$(OLIST) $(patsubst %.cc,%.o,$(CINTFILE))) -shared -o $(SHLIBFILE)

# Useful build targets
lib: $(LIBFILE) 
shlib: $(SHLIBFILE)
clean:
	rm -f $(OBJDIR)/*
	rm -f $(INCDIR)/*
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)

.PHONY : shlib lib default clean

-include $(addprefix $(INCDIR)/,$(DLIST))
