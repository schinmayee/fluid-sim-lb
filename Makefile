.PHONY: common projects
all: common projects

common:
	$(MAKE) -C common

projects:
	$(MAKE) -C projects

# remove object files
.PHONY: clean
clean:
	$(MAKE) -C common clean
	$(MAKE) -C projects clean
