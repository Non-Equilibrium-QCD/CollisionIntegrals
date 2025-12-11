# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fopenmp -g
LDFLAGS  := -lgsl -lgslcblas -lm -fopenmp -lfmt -lcuba

# Directories
SRC_DIR := src

# Output
TARGET := collint

# All source files
DEPS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/**/*.cpp)

# Default target
all: $(TARGET)

$(TARGET): $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/main.cpp $(LDFLAGS)

clean:
	rm -f $(TARGET) testint

testint: $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/testIntegral.cpp $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

debug: CXXFLAGS += -g -O0 -DDEBUG
debug: clean all

.PHONY: all clean run debug testint
