# Define the Fortran compiler and f2py command
FC = gfortran
F2PY = f2py3
F2PY_FLAGS = --f90flags="-fallow-argument-mismatch -std=legacy"
F2PY_FLAGS_77 = --f77flags="-fallow-argument-mismatch -std=legacy"
LIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# List of systems
SYSTEMS = ALC AlmHn Aln ArNO BrH2 BrHCl C2H2N2O C2H2O C2H6Cl C2H6F C2H6H C2H6O C2H6OH C2H7 C2N C3H7NO C6H5SH C6H6O C7H8S CH2O CH2O2 CH3NH2 CH3OH CH4Br \
         CH4Cl CH4CN CH4F CH4_H2O_H2O CH4O CH4OCl CH4OF CH4OH CH5 CH5O2 Cl2H ClH2 ClH2O ClNH3 COCO F2H F2H2 FH2 FH2O GeH4OH GeH5 H2CO H2O2 H2OBr H3 H3ClO H3p \
         H3S HClOH HO3 HOBr HOCO_Ar I2H IH2 K2Rb2 LiFH MCH MXH N2HOCp N2O N2O2 N3 N4 NaFH NaH2 NH3 NH3CN NH3OH NH4 NO2 O2H O2H3 O3 O4 OH2 OH3 SiH4Cl SiH5 SiOH YRH
# the existing systems that are not public for current version
# CH5p, CH4, H3Cl, H4O2

.PHONY: all clean check

all:
	@cd chempotpy; \
	for system in $(SYSTEMS); do \
		echo "Compiling system $$system..."; \
		cd $$system; \
		for surface in *.f90; do \
			if [ -f "$$surface" ]; then \
				name=$$(echo $$surface | sed 's/\.[^\.]*$$//'); \
				$(F2PY) -c -m $$name $(F2PY_FLAGS) $$surface $(LIBS) > $$name.log 2>&1; \
			fi; \
		done; \
                for surface in *.f; do \
                        if [ -f "$$surface" ]; then \
                                name=$$(echo $$surface | sed 's/\.[^\.]*$$//'); \
                                $(F2PY) -c -m $$name $(F2PY_FLAGS_77) $$surface $(LIBS) > $$name.log 2>&1; \
                        fi; \
                done; \
		cd ..; \
	done
	@echo "=========================================================================="
	@echo "=========================================================================="
	@echo "   "
	@echo "chempotpy compilation has been finished"
	@echo "   "
	@echo "=========================================================================="
	@echo "=========================================================================="

clean:
	@cd chempotpy; \
        for system in $(SYSTEMS); do \
                echo "Cleaning system $$system..."; \
                cd $$system; \
                rm -f *.so; \
                rm -f *.log; \
                rm -f *.o; \
                rm -f *.mod; \
                rm -f -r lib; \
                cd ..; \
        done


check:
	@cd chempotpy; \
	for system in $(SYSTEMS); do \
		echo "Checking system $$system..."; \
		cd $$system; \
		for surface in *.f *.f90; do \
			if [ -f "$$surface" ]; then \
				name=$$(echo $$surface | sed 's/\.[^\.]*$$//'); \
				success=$$(grep 'Removing build' $$name.log); \
				if [ ! -n "$$success" ]; then \
					echo "System $$system, surface $$surface did not compile successfully"; \
				fi; \
			fi; \
		done; \
		cd ..; \
	done
	@echo "=========================================================================="
	@echo "=========================================================================="
	@echo "   "
	@echo "Checking finished"
	@echo "   "
	@echo "=========================================================================="
	@echo "=========================================================================="

