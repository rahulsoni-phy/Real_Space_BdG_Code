# Platform: linux
EXENAME= SC_BdG_exe
### ------ Personal PC compilation ------------
CXX = g++
CPPFLAGS = -std=c++11
LDFLAGS = -llapack -lblas

CPPFLAGS += -Isrc
#CPPFLAGS += -g3
CPPFLAGS += -O3

#CPPFLAGS += -Isrc #-I/usr/include/eigen3
#CPPFLAGS += -I/usr/include/mkl/
LDFLAGS += # -lmkl_intel_lp64 -lmkl_core #-lmkl_intel_thread -lpthread -lm -ldl
LDFLAGS += #-fopenmp

STRIP_COMMAND = true

$(EXENAME): clean main.o
				$(CXX) $(CPPFLAGS) -o $(EXENAME)  main.o $(LDFLAGS)
				$(STRIP_COMMAND) $(EXENAME)

all: $(EXENAME)

clean:
				rm -f $(EXENAME) *.o

