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

