.PHONY: all clean

.DEFAULT: all


CXX = c++

CXXFLAGS = -ansi -pedantic -std=c++11  -O3

TARGETS = nqs_run

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

$(TARGETS) : %: main.cc
	$(CXX) $(CXXFLAGS) -o $@ $<
