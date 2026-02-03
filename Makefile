CC = g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

ASAN_FLAGS = -fsanitize=address -g -O1

CFLAGS = -Wall -std=c++17 $(ROOTCFLAGS) $(ASAN_FLAGS)
LDFLAGS = $(ROOTLIBS) $(ASAN_FLAGS)

SRC_DIR = src
INC_DIR = include
OBJ_DIR = build
BIN_DIR = bin

EXECUTABLE = $(BIN_DIR)/SABREsim

SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
DEPS = $(OBJECTS:.o=.d)

INCLUDES = -I$(INC_DIR)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -MMD -c $< -o $@

-include $(DEPS)
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)