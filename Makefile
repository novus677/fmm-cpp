CXX = clang++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -I./src
LDFLAGS = 

SOURCES = src/multipole.cpp src/local.cpp
OBJECTS = $(SOURCES:.cpp=.o)

NBODY_EXEC = test/nbody
SIMPLE_EXEC = test/simple

NBODY_OBJ = test/nbody.o
SIMPLE_OBJ = test/simple.o

ALL_OBJECTS = $(OBJECTS) $(NBODY_OBJ) $(SIMPLE_OBJ)
DEPS = $(ALL_OBJECTS:.o=.d)

all: nbody simple

nbody: $(NBODY_EXEC)
simple: $(SIMPLE_EXEC)

-include $(DEPS)

$(NBODY_EXEC): $(OBJECTS) $(NBODY_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(SIMPLE_EXEC): $(OBJECTS) $(SIMPLE_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c -o $@ $<

clean:
	rm -f $(ALL_OBJECTS) $(DEPS) $(NBODY_EXEC) $(SIMPLE_EXEC)

.PHONY: all clean nbody simple
