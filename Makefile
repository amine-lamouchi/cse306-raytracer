# Compiler
CXX = g++-13
# OpenMP Flag
OMPFLAG = -fopenmp
# Optimization flags for release mode
OPTFLAGS = -O3 -DNDEBUG
# Compiler flags
CXXFLAGS = $(OPTFLAGS) $(OMPFLAG) -std=c++17

# The build target executable
TARGET = main

# The source file(s)
SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	$(RM) $(TARGET)