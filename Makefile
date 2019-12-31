# Compiler stuff
COMPILER :=        g++-9

CFLAGS :=         -c -Wall -Wsign-compare -Wextra -fopenmp -O2 -std=c++11

INC_PATH :=       -iquote ./src \
                  -iquote ./test


# Object files for compiling + their paths
OBJS :=            Main.o \
					         Vector.o \
	                 Tensor.o \
									 Particle.o \
									 Body.o \
									 Neighbors.o \
									 Update.o \
									 Damage.o \
									 Contact.o \
									 Simulation.o \
									 Simulation_Setup.o \
									 VTK_File.o \
									 Data_Dump.o \
									 FEB_File.o

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
						 ./src/IO \
						 ./src/Simulation \
						 ./test




################################################################################
# Compile

compile: bin/SPH

bin/SPH: $(OBJ_PATHS)
	$(COMPILER) -fopenmp $(OBJ_PATHS) -o $@

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
obj/Particle.o: Particle.cc Particle.h Simulation.h Vector.h Tensor.h Errors.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Body class
obj/Body.o: Body.cc Body.h Simulation.h Vector.h Tensor.h Errors.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Neighbors.o: Neighbors.cc obj/Body.o
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Update.o: Update.cc obj/Body.o
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Damage.o: Damage.cc obj/Body.o
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Contact.o: Contact.cc obj/Body.o
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Simulation
obj/Simulation.o: Simulation.cc Simulation.h Body.h Particle.h Tensor.h Vector.h VTK_File.h Data_Dump.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Simulation_Setup.o: Simulation_Setup.cc Simulation.h Body.h  Particle.h Vector.h FEB_File.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# IO
obj/Data_Dump.o: Data_Dump.cc Data_Dump.h VTK_File.h Particle.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/VTK_File.o: VTK_File.cc VTK_File.h Data_Dump.h Particle.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/FEB_File.o: FEB_File.cc FEB_File.h Particle.h Body.h
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
