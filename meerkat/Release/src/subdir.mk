################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/bam_reader_test.cpp 

OBJS += \
./src/bam_reader_test.o 

CPP_DEPS += \
./src/bam_reader_test.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -static -pthread -std=c++11 -fopenmp -D__cplusplus=201103L -I/opt/gcc/4.8.5/include/c++/4.8.5 -I/opt/boost-1.57.0/include -I"/n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/src" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


