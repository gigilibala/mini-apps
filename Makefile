# COPYRIGHT @ Amin Hassani 2015

PROJECTS = stencil/1d benchmarks mantevo/minife-2.0/src mantevo/lulesh
LIBRARIES= ftlib comm-pat

TARGETS  = $(LIBRARIES) $(PROJECTS)
# TARGETS  = $(filter-out ./, $(dir $(shell find . -name 'Makefile' -print)))

BUILDDIRS = $(TARGETS:%=build-%)
INSTALLDIRS = $(TARGETS:%=install-%)
CLEANDIRS = $(TARGETS:%=clean-%)

all: $(BUILDDIRS)
$(TARGETS): $(BUILDDIRS)
$(BUILDDIRS):
	$(MAKE) -C $(@:build-%=%)

install: $(INSTALLDIRS)
$(INSTALLDIRS): 
	$(MAKE) -C $(@:install-%=%) install

clean: $(CLEANDIRS)
$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean


.PHONY: all clean install $(TARGETS) $(BUILDDIRS) $(INSTALLDIRS) $(CLEANDIRS) 
