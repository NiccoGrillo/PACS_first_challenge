# Define compiler variables
# CC = g++
# CFLAGS = -Wall -std=c++17
# LDFLAGS = -lgetpot

# # Specify output executable name
# OUTPUT = myprogram

# # List object (.o) files corresponding to each .cpp file
# SRCS_MAIN := $(wildcard main/*.cpp)
# SRCS_SRC := $(wildcard src/*.cpp)
# # OBJS_MAIN := $(addprefix obj/, $(notdir $(SRCS_MAIN)))
# # OBJS_SRC := $(addprefix obj/, $(notdir $(SRCS_SRC)))
# # OBJFILES := $(OBJS_MAIN) $(OBJS_SRC)
# OBJFILES = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) \
#           $(patsubst %.cpp,%.o,$(wildcard main/*.cpp)) 

# # Include paths for headers
# INCDIRS := ./include ./src
# INCPATHS := $(foreach dir, $(INCDIRS),-I$d)


# # Default target rule
# default: $(OUTPUT)

# # Linking stage
# $(OUTPUT): $(OBJFILES)
# 	$(CC) $(CFLAGS) $(INCPATHS) $^ -o $@ 

# # Compiling stage
# obj/%.o: main/%.cpp
# 	@mkdir -p obj
# 	@echo "Compiling $@"
# 	$(CC) $(CFLAGS) $(INCPATHS) -c $< -o $@

# obj/%.o: src/%.cpp
# 	@mkdir -p obj 
# 	@echo "Compiling $@"
# 	$(CC) $(CFLAGS) $(INCPATHS) -c $< -o $@

# # Phony targets
# .PHONY: clean
# clean:
# 	rm -f $(OBJFILES) $(OUTPUT)

# ifneq ($(MAKECMDGOALS),clean)
# -include $(shell find dep -name '*.d')
#endif

















############        OLD_VERSION        ###########




CC = g++
CFLAGS = -Wall -std=c++17
# LDFLAGS = -lgetpot

# Specify output executable name
OUTPUT = myprogram

# List object (.o) files corresponding to each .cpp file
CPP_MAIN = $(wildcard main/*.cpp)
CPP_SRC = $(wildcard src/*.cpp)
OBJECTS_MAIN = $(patsubst %.cpp,%.o,$(wildcard main/*.cpp)) 
OBJECTS_SRC =  $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) # Add this line
# Include paths for headers
INCLUDEPATHS = -I./include -I./src

# Default target rule
$(OUTPUT): $(OBJECTS_MAIN)
	@echo "yoyoyoyo"
	$(CC) $(CFLAGS) $^ -o $@ 

# Compile objects from sources
src/%.o : src/%.cpp
	@echo "compiling $@ from files $<"
	$(CC) $(CFLAGS) $(INCLUDEPATHS) -c $< -o $@
main/%.o : $(OBJECTS_SRC)
	@echo "compiling $@ from files $<"
	$(CC) $(CFLAGS) $(INCLUDEPATHS)  $< -c $(CPP_MAIN) -o $@

clean:
	rm -f $(OBJECTS_MAIN) $(OBJECTS_SRC) $(OUTPUT)