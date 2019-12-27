# Compiler stuff
COMPILER :=        g++-9

CFLAGS :=         -c -Wall -Wsign-compare -Wextra -O2 -std=c++11

INC_PATH :=       -iquote ./src \
                  -iquote ./test


# Object files for compiling + their paths
OBJS :=            Main.o \
					         Vector.o \
	                 Tensor.o \
									 Particle.o \
									 Particle_Contact.o \
									 Particle_Damage.o \
									 Particle_Debugger.o \
									 Particle_Neighbors.o \
									 Particle_Update.o \
									 Body.o \
									 Simulation.o

OBJ_PATHS :=       $(patsubst %,obj/%,$(OBJS))


# Object files for tests + their paths
TEST_OBJS :=       Vector.o \
                   Tensor.o \
									 Tests.o
TEST_OBJ_PATHS := $(patsubst %,obj/%,$(TEST_OBJS))



# Tells Make where to search for dependencies
VPATH :=     ./bin ./obj \
	           ./src \
             ./src/Vector ./src/Tensor \
						 ./src/Particle ./src/Body \
						 ./src/Simulation \
						 ./test




################################################################################
# Compile

compile: bin/SPH

bin/SPH: $(OBJ_PATHS)
	$(COMPILER) $(OBJ_PATHS) -o $@

obj/Main.o: Main.cc
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@




################################################################################
# Test

test: bin/Test

bin/Test: $(TEST_OBJ_PATHS)
	$(COMPILER) $(TEST_OBJ_PATHS) -o $@

obj/Tests.o: Tests.cc Vector_Tests.cc Tensor_Tests.cc
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@




################################################################################
# Object files

# Vector class.
obj/Vector.o: Vector.cc Vector.h Tensor.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Tensor class.
obj/Tensor.o: Tensor.cc Tensor.h Vector.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Particle Class
obj/Particle.o: Particle.cc Particle.h Vector.h Tensor.h Errors.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Particle_Update.o: Particle_Update.cc Particle.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Particle_Damage.o: Particle_Damage.cc Particle.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Particle_Contact.o: Particle_Contact.cc Particle.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Particle_Debugger.o: Particle_Debugger.cc Particle.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Body class
obj/Body.o: Body.cc Body.h Vector.h Tensor.h Errors.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Neighbors.o: Neighbors.cc Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



obj/Simulation.o: Simulation.cc Simulation.h Vector.h Tensor.h Particle.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



obj/Data_Dump.o: Data_Dump.cc Data_Dump.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



################################################################################
# Clean and help

clean:
	rm ./obj/*.o
	rm ./bin/SPH


help:
	$(info )
	$(info Commands for SPH Makefile:)
	$(info make           - compiles code to 'SPH' (in /bin))
	$(info make compile   - same as 'make')
	$(info make test      - makes unit test file (in /bin))
	$(info make clean     - clears all object files (in /obj) and binary files (in \bin))
	$(info make help      - displays this message)
	$(info )
