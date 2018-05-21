CXX=g++
LD=g++

LDFLAGS := $(shell root-config --glibs )
CFLAGS := -O3 $(shell root-config --cflags )

SRC_FILES := $(wildcard *.cxx)
EXE_FILES := $(SRC_FILES:%.cxx=%)

.PHONY: all clean

all: $(EXE_FILES)

%: %.o
	$(LD) $(LDFLAGS) $^ -o $@

%.o: %.cxx
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXE_FILES:%=%.o) $(EXE_FILES)
