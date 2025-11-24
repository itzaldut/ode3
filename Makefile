# top level Makefile
# this "simply" runs the make in the two subdirectories: odelib and src


# top level Makefile
# this "simply" runs the make in the two subdirectories: odelib and src


ODELIB =  $(PWD)/odelib
export ODELIB

SUBDIRS = odelib src

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS): 
	make -C $@

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean ;\
	done
