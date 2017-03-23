################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/meerkat/cluster/AlternativeMapper.cpp \
../src/meerkat/cluster/Cluster.cpp \
../src/meerkat/cluster/ClusterContainer.cpp \
../src/meerkat/cluster/ClusterSelector.cpp \
../src/meerkat/cluster/Read.cpp \
../src/meerkat/cluster/ReadInfo.cpp \
../src/meerkat/cluster/Readpair.cpp 

OBJS += \
./src/meerkat/cluster/AlternativeMapper.o \
./src/meerkat/cluster/Cluster.o \
./src/meerkat/cluster/ClusterContainer.o \
./src/meerkat/cluster/ClusterSelector.o \
./src/meerkat/cluster/Read.o \
./src/meerkat/cluster/ReadInfo.o \
./src/meerkat/cluster/Readpair.o 

CPP_DEPS += \
./src/meerkat/cluster/AlternativeMapper.d \
./src/meerkat/cluster/Cluster.d \
./src/meerkat/cluster/ClusterContainer.d \
./src/meerkat/cluster/ClusterSelector.d \
./src/meerkat/cluster/Read.d \
./src/meerkat/cluster/ReadInfo.d \
./src/meerkat/cluster/Readpair.d 


# Each subdirectory must supply rules for building sources it contributes
src/meerkat/cluster/%.o: ../src/meerkat/cluster/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -static -pthread -std=c++11 -fopenmp -D__cplusplus=201103L -I/opt/gcc/4.8.5/include/c++/4.8.5 -I/opt/boost-1.57.0/include -I"/n/data1/bch/genetics/lee/kyu/eclipse/workspaces/genomics/mybamtools/src" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


