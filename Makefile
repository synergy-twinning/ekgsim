.PHONY: clean All

All:
	@echo "----------Building project:[ ekgSim - Release ]----------"
	@"$(MAKE)" -f  "ekgSim.mk"
clean:
	@echo "----------Cleaning project:[ ekgSim - Release ]----------"
	@"$(MAKE)" -f  "ekgSim.mk" clean
