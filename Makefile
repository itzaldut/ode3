# top level Makefile
# this "simply" runs the make in the two subdirectories: odelib and src


# Compiler and flags
GXX = g++
CXXFLAGS = -std=c++17 -Wall
LDFLAGS =

# Paths
SRC_DIR = src
ODELIB = odelib

# Files
TARGET = baseball1
SRC = $(SRC_DIR)/$(TARGET).cpp
HEADERS = $(ODELIB)/RKn.hpp

# Default target
all: $(TARGET)

# Build rule
$(TARGET): $(SRC) $(HEADERS)
	$(GXX) $(CXXFLAGS) -I$(ODELIB) -o $@ $(SRC) $(LDFLAGS)

# Clean rule
clean:
	rm -f $(TARGET)

