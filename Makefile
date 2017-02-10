.PHONY: clean All

All:
	@echo "----------Building project:[ simlib - Debug ]----------"
	@cd "simlib" && "$(MAKE)" -f  "simlib.mk"
	@echo "----------Building project:[ ekgSim - Debug ]----------"
	@"$(MAKE)" -f  "ekgSim.mk"
clean:
	@echo "----------Cleaning project:[ simlib - Debug ]----------"
	@cd "simlib" && "$(MAKE)" -f  "simlib.mk"  clean
	@echo "----------Cleaning project:[ ekgSim - Debug ]----------"
	@"$(MAKE)" -f  "ekgSim.mk" clean
