##########################################################################
# Build the tar.gz for building with configure/make.                     #
##########################################################################

NAME ?= PhreeqcRM

UPDATE := $(shell svn up ..)

REL := $(shell svn st -v ../configure.ac | cut -b 10- | ( read one two ; echo $$one ) )

TARBALL := $(NAME)-$(VER)-$(REL).tar.gz

ZIPBALL := $(NAME)-$(VER)-$(REL).zip

REPOS_PATH ?= trunk

all : $(TARBALL) $(ZIPBALL)

$(TARBALL) : dist.sh
ifeq ($(DATE),)
	$(error usage: make DATE=<m/d/yyyy> VER=<x.y.z> [REL=x])
endif
ifeq ($(VER),)
	$(error usage: make DATE=<m/d/yyyy> VER=<x.y.z> [REL=x])
endif
	@echo "Creating $(TARBALL)"
	@rm -f $(TARBALL)
	time ./dist.sh -v $(VER) -r $(REL) -d $(DATE) -pr $(REPOS_PATH) 2>&1 | tee dist-$(REL).log


$(ZIPBALL) : $(TARBALL)
	@echo "Creating $(ZIPBALL)"
	@rm -f $(ZIPBALL)
	time ./dist.sh -v $(VER) -r $(REL) -d $(DATE) -pr $(REPOS_PATH) -win 2>&1 | tee distw-$(REL).log


usage:
	@echo ""
	@echo "usage: make DATE=<m/d/yyyy> VER=<x.y.z> [REL=x]"
	@echo "    ie make DATE=8/12/2010 VER=2.17.3"
	@echo "    ie make DATE=8/12/2010 VER=2.17.3 REL=4000"
	@echo ""
