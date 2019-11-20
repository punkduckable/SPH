#Variables
COMPILER := g++-9

CFLAGS := -c -Wall -Wsign-compare -Wextra -O2 -std=c++11

INC_PATH :=   -iquote ./src \
              -iquote ./test

OBJS :=        Main.o \
					     Vector.o \
	             Tensor.o

PATH_OBJS := $(patsubst %,obj/%,$(OBJS))

VPATH :=     ./bin ./obj ./src \
             ./src/Vector ./src/Tensor \
						 ./test

# All
all: $(PATH_OBJS) bin/SPH



# General rules
bin/SPH: $(PATH_OBJS)
	$(COMPILER) $(PATH_OBJS) -o $@

obj/Main.o: Main.cc
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Rules for the Vector class.
obj/Vector.o: Vector.cc Vector.h Tensor.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Rules for the Tensor class.
obj/Tensor.o: Tensor.cc Tensor.h Vector.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



clean:
	rm ./obj/*.o
	rm ./bin/SPH


help:
	$(info )
	$(info Commands for SPH Makefile:)
	$(info make       - compiles code to 'SPH' (in \bin))
	$(info make all   - same as 'make')
	$(info make clean - clears all object files (in \obj) and binary files (in \bin))
	$(info make help  - displays this message)
	$(info )
