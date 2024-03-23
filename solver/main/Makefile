CPPFLAGS=-I../src -I../include
CXXFLAGS=-std=c++17 -Wall -O2
#OBJECTS = Utilities.o optimizer.o

all: main
main: main.o
	g++ main.o -o main

main.o: solver.o
	g++ -l solver.o -c main.cpp $(CPPFLAGS)
solver.o: ../src/solver.cpp
	g++ -c ../src/solver.cpp $(CPPFLAGS)
clean:
	rm main main.o solver.o











# # CC = g++
# # CFLAGS = -Wall -std=c++17

# # # Specify output executable name
# # OUTPUT = myprogram

# # # List object (.o) files corresponding to each .cpp file
# # OBJECTS = $(patsubst %.cpp,%.o,$(wildcard *.cpp)) \
# # 			$(patsubst %.cpp,%.o,$(wildcard ../src/*.cpp))
# #            # Add this line
# # # Include paths for headers
# # INCLUDEPATHS = -I../include -I../src

# # # Default target rule
# # $(OUTPUT): $(OBJECTS)
# # 	$(CC) $(CFLAGS) $^ -o $@ 

# # # Compile objects from sources
# # ../src/%.o : ../src/%.cpp
# # 	@echo "compiling $@ from files $<"
# # 	$(CC) $(CFLAGS) $(INCLUDEPATHS) -c $< -o $@
# # %.o : %.cpp 
# # 	@echo "compiling $@ from files $<"
# # 	$(CC) $(CFLAGS) $(INCLUDEPATHS) -c $< -o $@

# # clean:
# # 	rm -f $(OBJECTS) $(OUTPUT)



# CC = g++
# CFLAGS = -Wall -std=c++17 -MMD -MP

# # Include paths for headers
# INCLUDEPATHS = -I../include -I../src

# # Specify output executable name
# OUTPUT = myprogram

# # Source directories
# SRCDIR = . ../src

# # Find all .cpp files and their corresponding .o files
# SOURCES := $(foreach dir,$(SRCDIR),$(wildcard $(dir)/*.cpp))
# OBJECTS := $(SOURCES:.cpp=.o)
# DEPS := $(patsubst %.cpp,%.o,$(wildcard ../src/*.cpp))

# # Default target rule
# $(OUTPUT): $(OBJECTS)
# 	$(CC) $(CFLAGS) $^ -o $@ $(INCLUDEPATHS)

# # Compile objects from sources
# %.o: %.cpp
# 	@echo "Compiling $< to $@"
# 	$(CC) $(CFLAGS) $(INCLUDEPATHS) -c $< -o $@

# -include $(DEPS)

# clean:
# 	rm -f $(OBJECTS) $(DEPS) $(OUTPUT)
