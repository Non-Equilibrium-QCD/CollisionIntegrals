# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -fopenmp -g -march=native
CXXFLAGS += -O3 -funroll-loops -flto
CXXFLAGS += -I./third_party/fmt/include/ -L./third_party/fmt/build -lfmt
LDFLAGS  := -lgsl -lgslcblas -lm -lfmt -fopenmp -lcuba

BOLDGREEN := \033[1;32m
RESET := \033[0m

# Directories
SRC_DIR := src
TEST_DIR := test

# Output
TARGET := collint

# All source files
DEPS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/**/*.cpp)
TEST_DEPS := $(wildcard $(TEST_DIR)/*.cpp) $(DEPS)

# Default target
all: $(TARGET)

$(TARGET): $(DEPS) libfmt
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/main.cpp $(LDFLAGS)

clean:
	rm -f $(TARGET) testint testsimple

testint: $(TEST_DEPS) libfmt
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_DIR)/testIntegral.cpp $(LDFLAGS)

testsimple: $(TEST_DEPS) libfmt
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_DIR)/testSimpleIntegral.cpp $(LDFLAGS)
	./testsimple

test: testint testsimple
	./testint
	./testsimple

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

# Examples: make example NAME=YM
EXAMPLES_DIR := examples
AVAILABLE_EXAMPLES := $(basename $(notdir $(wildcard $(EXAMPLES_DIR)/*.cpp)))

example-help:
	@echo "Usage: make example NAME=<name>"
	@echo "       make example-run NAME=<name>"
	@echo ""
	@echo "Available examples:"
	@for ex in $(AVAILABLE_EXAMPLES); do echo "  - $$ex"; done
	@echo ""
	@echo "Example: make example NAME=YM"

example:
ifndef NAME
	@$(MAKE) --no-print-directory example-help
	@exit 1
else ifeq ($(wildcard $(EXAMPLES_DIR)/$(NAME).cpp),)
	@echo "Error: $(EXAMPLES_DIR)/$(NAME).cpp not found"
	@echo ""
	@$(MAKE) --no-print-directory example-help
	@exit 1
else
	@echo -e "${BOLDGREEN}Building example $(NAME) successfully.${RESET}"
	$(CXX) $(CXXFLAGS) -o $(EXAMPLES_DIR)/$(NAME) $(EXAMPLES_DIR)/$(NAME).cpp $(LDFLAGS)
endif

example-run: example libfmt
	@echo -e "${BOLDGREEN}Running example $(NAME)...${RESET}"
	./$(EXAMPLES_DIR)/$(NAME)

libfmt: third_party/fmt/build/libfmt.a

third_party/fmt/build/libfmt.a:
	cd third_party/fmt && mkdir -p build && cd build &&\
		cmake -DFMT_DOC=OFF -DFMT_TEST=OFF -DFMT_INSTALL=OFF -DFMT_MASTER_PROJECT=OFF .. &&\
		make -j

.PHONY: all clean run debug testint testsimple test example example-help example-run docs docs-xml docs-pdf docs-open docs-pdf-open docs-clean docs-all docs-watch
