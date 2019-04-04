.PHONY: all clean

CFLAGS += -std=c++11
# CFLAGS += -g
CFLAGS += -O3
LDFLAGS += -lm     # link to math library

all: train test

train: train.cpp hmm.h
	g++ -o $@ train.cpp $(CFLAGS)

test: test.cpp hmm.h
	g++ -o $@ test.cpp $(CFLAGS)

clean:
	$(RM) train test
