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
									 Body.o Neighbors.o Update.o Damage.o Contact.o Body_IO.o \
									 Simulation.o Simulation_Setup.o Timing.o Boundary_Conditions.o \
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
	$(COMPILER) -fopenmp $(OBJ_PATHS) -o $@

obj/Main.o: Main.cc
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@




################################################################################
# Test

test: bin/Test

bin/Test: $(TEST_OBJ_PATHS)
	$(COMPILER) $(TEST_OBJ_PATHS) -o $@

obj/Tests.o: Tests.cc Vector_Tests.cc Tensor_Tests.cc IO_Tests.cc
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

obj/Simulation_Setup.o: Simulation_Setup.cc Simulation.h Body.h Particle.h Vector.h FEB_File.h Errors.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Timing.o: Timing.cc Simulation.h
	$(COMPILER) $(CFLAGS) $(INC_PATH) $< -o $@

obj/Boundary_Conditions.o: Boundary_Conditions.cc Simulation.h Body.h Particle.h Vector.h Array.h
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
