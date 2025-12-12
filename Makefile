# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fopenmp -g -march=native
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

# Documentation
docs:
	doxygen Doxyfile
	@echo "Documentation generated at docs/html/index.html"

docs-pdf: docs
	$(MAKE) -C docs/latex || true
	@cp docs/latex/refman.pdf docs/CollIntegral.pdf 2>/dev/null || { echo "PDF generation failed. Install texlive: sudo apt install texlive-latex-extra"; exit 1; }
	@echo "PDF generated at docs/CollIntegral.pdf"

docs-open: docs
	xdg-open docs/html/index.html 2>/dev/null || open docs/html/index.html

docs-pdf-open: docs-pdf
	xdg-open docs/CollIntegral.pdf 2>/dev/null || open docs/CollIntegral.pdf

docs-clean:
	rm -rf docs/html docs/latex docs/xml docs/CollIntegral.pdf

.PHONY: all clean run debug testint docs docs-pdf docs-open docs-pdf-open docs-clean
