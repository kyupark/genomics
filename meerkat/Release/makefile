################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/third/subdir.mk
-include src/meerkat/srd/subdir.mk
-include src/meerkat/preprocess/subdir.mk
-include src/meerkat/mpd/subdir.mk
-include src/meerkat/mechanism/subdir.mk
-include src/meerkat/filter/subdir.mk
-include src/meerkat/dre/subdir.mk
-include src/meerkat/cns/subdir.mk
-include src/meerkat/cluster/subdir.mk
-include src/meerkat/blast_bp/subdir.mk
-include src/meerkat/alg/subdir.mk
-include src/meerkat/subdir.mk
-include src/castle/subdir.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: meerkat

dependents:
	-cd /n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/Release && $(MAKE) clean all

# Tool invocations
meerkat: $(OBJS) $(USER_OBJS) /n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/Release/mybamtools
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++ -static -pthread -std=c++11 -fopenmp -L"/n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/lib" -L"/n/data1/hms/dbmi/park/alee/kyu/eclipse/workspaces/genomics/lib" -L/opt/boost-1.57.0/lib/ -o "meerkat" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS) meerkat
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
/n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/Release/mybamtools:

-include ../makefile.targets
