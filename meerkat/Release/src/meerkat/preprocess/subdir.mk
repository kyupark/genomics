################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/meerkat/preprocess/Histogram.cpp \
../src/meerkat/preprocess/ParallelBamReader.cpp \
../src/meerkat/preprocess/ReadGroup.cpp \
../src/meerkat/preprocess/ReadGroupAlt.cpp 

OBJS += \
./src/meerkat/preprocess/Histogram.o \
./src/meerkat/preprocess/ParallelBamReader.o \
./src/meerkat/preprocess/ReadGroup.o \
./src/meerkat/preprocess/ReadGroupAlt.o 

CPP_DEPS += \
./src/meerkat/preprocess/Histogram.d \
./src/meerkat/preprocess/ParallelBamReader.d \
./src/meerkat/preprocess/ReadGroup.d \
./src/meerkat/preprocess/ReadGroupAlt.d 


# Each subdirectory must supply rules for building sources it contributes
src/meerkat/preprocess/%.o: ../src/meerkat/preprocess/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -static -pthread -std=c++11 -fopenmp -D__cplusplus=201103L -I/opt/gcc/4.8.5/include/c++/4.8.5 -I/opt/boost-1.57.0/include -I"/n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/src" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


