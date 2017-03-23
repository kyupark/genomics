################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/third/faidx.c \
../src/third/razf.c 

CPP_SRCS += \
../src/third/gzstream.cpp \
../src/third/xfaidx.cpp 

OBJS += \
./src/third/faidx.o \
./src/third/gzstream.o \
./src/third/razf.o \
./src/third/xfaidx.o 

C_DEPS += \
./src/third/faidx.d \
./src/third/razf.d 

CPP_DEPS += \
./src/third/gzstream.d \
./src/third/xfaidx.d 


# Each subdirectory must supply rules for building sources it contributes
src/third/%.o: ../src/third/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -static -pthread -fopenmp -I/n/data1/hms/dbmi/park/vincent_lim/codes/Meerkat/src/mybamtools/src -I/opt/boost-1.57.0/include -I/opt/gcc/4.8.5/include/c++/4.8.5 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/third/%.o: ../src/third/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -static -pthread -std=c++11 -fopenmp -D__cplusplus=201103L -I/opt/gcc/4.8.5/include/c++/4.8.5 -I/opt/boost-1.57.0/include -I"/n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/src" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


