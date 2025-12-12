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

# Documentation (Sphinx + Breathe + Doxygen)
SPHINXBUILD   = sphinx-build
SPHINXSRCDIR  = docs/source
SPHINXBLDDIR  = docs/build

docs-xml:
	cd docs && doxygen Doxyfile

docs: docs-xml
	$(SPHINXBUILD) -b html $(SPHINXSRCDIR) $(SPHINXBLDDIR)/html
	@echo "Documentation generated at $(SPHINXBLDDIR)/html/index.html"

docs-pdf: docs-xml
	$(SPHINXBUILD) -b latex $(SPHINXSRCDIR) $(SPHINXBLDDIR)/latex
	$(MAKE) -C $(SPHINXBLDDIR)/latex all-pdf || true
	@cp $(SPHINXBLDDIR)/latex/CollIntegral.pdf docs/ 2>/dev/null || { echo "PDF generation failed. Install texlive."; exit 1; }
	@echo "PDF generated at docs/CollIntegral.pdf"

docs-open: docs
	xdg-open $(SPHINXBLDDIR)/html/index.html 2>/dev/null || open $(SPHINXBLDDIR)/html/index.html

docs-pdf-open: docs-pdf
	xdg-open docs/CollIntegral.pdf 2>/dev/null || open docs/CollIntegral.pdf

docs-clean:
	rm -rf $(SPHINXBLDDIR) docs/xml docs/CollIntegral.pdf

docs-all: docs docs-pdf
	@echo "HTML and PDF documentation generated"

docs-watch:
	@echo "Watching docs/source for changes... (Ctrl+C to stop)"
	@while true; do \
		inotifywait -qre modify,create,delete $(SPHINXSRCDIR) 2>/dev/null || \
		fswatch -1 $(SPHINXSRCDIR) 2>/dev/null || \
		{ echo "Install inotify-tools or fswatch: sudo pacman -S inotify-tools"; exit 1; }; \
		echo "Change detected, rebuilding..."; \
		$(MAKE) docs-all; \
		echo "Done. Waiting for changes..."; \
	done

.PHONY: all clean run debug testint docs docs-xml docs-pdf docs-open docs-pdf-open docs-clean docs-all docs-watch
