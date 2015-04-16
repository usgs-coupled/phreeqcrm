##########################################################################
# Build the tar.gz for building with configure/make.                     #
##########################################################################

ifeq ($(DATE),)
    $(error usage: make DATE=<m/d/yyyy> VER=<x.y.z> [REL=x])
endif
ifeq ($(VER),)
    $(error usage: make DATE=<m/d/yyyy> VER=<x.y.z> [REL=x])
endif


NAME ?= PhreeqcRM

LNAME := $(shell echo $(NAME) | tr '[A-Z]' '[a-z]')

UPDATE := $(shell svn up ../..)

REL := $(shell svn st -v ../../configure.ac | cut -b 10- | ( read one two ; echo $$one ) )

CONF_PKG := $(LNAME)-$(VER)-$(REL).tar.gz

CMAKE_PKG := $(LNAME)-$(VER)-$(REL).zip

TARBALL  := $(NAME)-$(VER)-$(REL).tar.gz

ZIPBALL  := $(NAME)-$(VER)-$(REL).zip

REPOS_PATH ?= trunk



all : $(CONF_PKG) $(CMAKE_PKG)
	@echo "DATE=$(DATE)"        >> history
	@echo "VER=$(VER)"          >> history
	@echo "REL=$(REL)"          >> history
	@echo "===================" >> history


$(CONF_PKG) : ../$(TARBALL)
	rm -rf $(NAME)-$(VER)-$(REL)
	tar xvzf ../$(TARBALL)
	cd $(NAME)-$(VER)-$(REL) && autoreconf -fvi
	cd $(NAME)-$(VER)-$(REL) && mkdir _build
	cd $(NAME)-$(VER)-$(REL) && cd _build && ../configure
	cd $(NAME)-$(VER)-$(REL) && cd _build && $(MAKE) distcheck
	mv $(NAME)-$(VER)-$(REL)/_build/$(CONF_PKG) .
	rm -rf $(NAME)-$(VER)-$(REL)

$(CMAKE_PKG) : ../$(ZIPBALL) $(CONF_PKG)
	rm -rf $(NAME)-$(VER)-$(REL)
	unzip ../$(ZIPBALL)
	cd $(NAME)-$(VER)-$(REL) && autoreconf -fvi
	cd $(NAME)-$(VER)-$(REL) && mkdir _build
	cd $(NAME)-$(VER)-$(REL) && cd _build && ../configure
	cd $(NAME)-$(VER)-$(REL) && cd _build && $(MAKE) dist-zip
	mv $(NAME)-$(VER)-$(REL)/_build/$(CMAKE_PKG) .
	rm -rf $(NAME)-$(VER)-$(REL)

usage:
	@echo ""
	@echo "usage: make DATE=<m/d/yyyy> VER=<x.y.z> [REL=x]"
	@echo "    ie make DATE=8/12/2010 VER=2.17.3"
	@echo "    ie make DATE=8/12/2010 VER=2.17.3 REL=4000"
	@echo ""
