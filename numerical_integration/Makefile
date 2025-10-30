EIGEN_INC ?= $(abspath eigen-3.4.0)
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -I$(EIGEN_INC)
LDFLAGS ?= -lm

SRC = basic_toolbox.cpp \
      generate_graph_from_sequence_toolbox.cpp \
      ecosystem_numerical_integration_toolbox.cpp

DRIVERS = Ecosystem_Dynamics_on_Graphs_in_Temperature \
          Ecosystem_Dynamics_on_RRG_zero_Temperature \
          Ecosystem_Dynamics_on_DirectedER_zero_Temperature

BIN_DIR = ../build/bin

all: $(DRIVERS:%=$(BIN_DIR)/%)

$(BIN_DIR)/%: $(SRC) $(@F).cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SRC) $(@F).cpp -o $@ $(LDFLAGS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)

.PHONY: all clean
