CXX = g++
CXXFLAGS = -std=c++17 -O2
EIGEN ?= /usr/include/eigen3

SRCS = main.cc Molecule.cc Hessian.cc Mass.cc VibrationalAnalysis.cc
OBJS = $(SRCS:.cc=.o)
TARGET = vib

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -I$(EIGEN) -o $(TARGET)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -I$(EIGEN) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
