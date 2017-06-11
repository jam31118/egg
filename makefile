# Required Program for Compilation
CC = gcc

# Compile Configuration
CFLAGS = -O2 -Wall --std=c11 -I./include -lm

# Required Source Code for Compilation
COMOBJ = ./lib/config.c 

# Resulting Executive Binary files by Compilation
TARGET1 = even
TARGET2 = prop

# Resulting Data files after execution of the binaries
DATA = energy.dat tdse.dat y.dat psi0.dat

# Define what to do when "make" without arguments
all: $(TARGET1) $(TARGET2)

# Define for each arguments for "make"
$(TARGET1): $(COMOBJ) ./lib/$(TARGET1).c
	$(CC) -o $@ $^ $(CFLAGS)

prop: $(COMOBJ) ./lib/$(TARGET2).c
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	$(RM) $(TARGET1) $(TARGET2) $(DATA)
