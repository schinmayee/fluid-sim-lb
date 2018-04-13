.PHONY: flip-current-lb flip-lb particles-to-levelset
all: flip-current-lb flip-lb particles-to-levelset

flip-current-lb:
	$(MAKE) -C flip-current-lb

flip-lb:
	$(MAKE) -C flip-lb

particles-to-levelset:
	$(MAKE) -C particles-to-levelset

# remove object files
.PHONY: clean
clean:
	$(MAKE) -C flip-current-lb clean
	$(MAKE) -C flip-lb clean
	$(MAKE) -C particles-to-levelset clean
