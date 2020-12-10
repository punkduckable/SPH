# Compiler stuff
COMPILER :=        g++-10

CFLAGS   :=        -c -Wall -Wsign-compare -Wextra -fopenmp -O2 -std=c++11 -DNDEBUG

LNFLAGS  :=

INC_PATH :=        -iquote ./src \
                   -iquote ./test


# Object files for compiling + their paths
OBJS :=            Main.o \
					         Vector.o \
	                 Tensor.o \
									 Particle.o \
									 Body.o Neighbors.o Update.o Damage.o Contact.o Body_IO.o \
									 Simulation.o Simulation_Setup.o Timing.o Boundary_Conditions.o Setup_File_Reader.o \
									 Save_Simulation.o Load_Simulation.o FEB_File.o IO_Ops.o

OBJ_PATHS :=       $(patsubst %,obj/%,$(OBJS))


# Object files for tests + their paths
TEST_OBJS :=       Vector.o \
                   Tensor.o \
									 IO_Ops.o \
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
	$(COMPILER) $(LNFLAGS) -fopenmp $(OBJ_PATHS) -o $@



################################################################################
# Test

test: bin/Test

bin/Test: $(TEST_OBJ_PATHS)
	$(COMPILER) $(TEST_OBJ_PATHS) -o $@



################################################################################
# Debug

debug: CFLAGS := -c -Wall -Wsign-compare -Wextra -g -Og -fno-inline -std=c++11
debug: bin/Debug

bin/Debug: $(OBJ_PATHS)
		$(COMPILER) $(OBJ_PATHS) -o $@



################################################################################
# Profiling

prof: CFLAGS := -c -Wall -Wsign-compare -Wextra -pg -O0 -fno-inline -std=c++11 -fopenmp
prof: LNFLAGS := -pg
prof: bin/Prof

bin/Prof: $(OBJ_PATHS)
		$(COMPILER) $(LNFLAGS) -fopenmp $(OBJ_PATHS) -o $@



################################################################################
# Object files

# Main
obj/Main.o: Main.cc
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Test
obj/Tests.o: Tests.cc Vector_Tests.cc Tensor_Tests.cc IO_Tests.cc List_Tests.cc Array_Tests.cc Body_Tests.cc List.h Array.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Vector class.
obj/Vector.o: Vector.cc Vector.h Tensor.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Tensor class.
obj/Tensor.o: Tensor.cc Tensor.h Vector.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Particle Class
obj/Particle.o: Particle.cc Particle.h Simulation.h Vector.h Tensor.h Errors.h Body.h Load_Simulation.h Save_Simulation.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Body class
obj/Body.o: Body.cc Body.h Vector.h Particle.h Load_Simulation.h Save_Simulation.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Neighbors.o: Neighbors.cc Body.h Particle.h Vector.h List.h Array.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Update.o: Update.cc Body.h Simulation.h Particle.h Vector.h List.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Damage.o: Damage.cc Body.h Particle.h Vector.h List.h Array.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Contact.o: Contact.cc Body.h Simulation.h Particle.h Vector.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Body_IO.o: Body_IO.cc Body.h Particle.h Vector.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# Simulation
obj/Simulation.o: Simulation.cc Simulation.h Body.h Particle.h Vector.h Load_Simulation.h Save_Simulation.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Simulation_Setup.o: Simulation_Setup.cc Simulation.h Body.h Particle.h Vector.h FEB_File.h Array.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Timing.o: Timing.cc Simulation.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Boundary_Conditions.o: Boundary_Conditions.cc Simulation.h Body.h Particle.h Vector.h Array.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Setup_File_Reader.o: Setup_File_Reader.cc Simulation.h IO_Ops.h Body.h Vector.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



# IO
obj/Save_Simulation.o: Save_Simulation.cc Save_Simulation.h Body.h Particle.h Tensor.h Vector.h IO_Ops.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Load_Simulation.o: Load_Simulation.cc Load_Simulation.h Body.h Particle.h Tensor.h Vector.h IO_Ops.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/FEB_File.o: FEB_File.cc FEB_File.h Particle.h Body.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/IO_Ops.o: IO_Ops.cc IO_Ops.h Vector.h Simulation.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@



################################################################################
# Clean and help

clean:
	rm ./obj/*.o
	rm ./bin/*



help:
	$(info )
	$(info Commands for SPH Makefile:)
	$(info make           - compiles code to 'SPH' (in /bin))
	$(info make compile   - same as 'make')
	$(info make test      - makes unit test file (in /bin))
	$(info make debug     - makes Debug file for debugging (in /bin))
	$(info make clean     - clears all object files (in /obj) and binary files (in \bin))
	$(info make help      - displays this message)
	$(info )
