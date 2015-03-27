OPT_OR_DEBUG = -O3

CXX = g++
CXXFLAGS := -Wall -Wshadow $(OPT_OR_DEBUG)
LINKFLAGS := 

OFILES=\
	experiment.o

PROG = experiment

all: $(PROG)

$(PROG): $(OFILES)
	$(CXX) $(CXXFLAGS) $(OFILES) -o $@ $(LINKFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -MMD -MP -MF $(<:%.cc=.%.d) -c $<

.PHONY: clean
clean:
	rm -f *~ *.o $(PROG)

.PHONY: cleanall
cleanall: clean
	rm -f .*.d

# Import dependences
-include $(OFILES:%.o=.%.d)
-include .unittests.d
