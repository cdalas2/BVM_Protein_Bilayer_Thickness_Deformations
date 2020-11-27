#website that taught me to make this make file: http://nuclear.mutantstargoat.com/articles/make/
src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
dep = $(obj:.o=.d)  # one dependency file for each source
prefix = clover5
suffix = 030
conda_path = /home/calas/anaconda3/envs/arblib_env/

CXXFLAGS = -MMD -fopenmp -O3 -I$(conda_path)include -I$(PWD) -L$(conda_path)lib# option to generate a .d file during compilation
LDFLAGS = -larb -lflint -lmpfr -lgmp -lntl -lgf2x
CC = gcc
C++ = g++
$(prefix)eps$(suffix): $(obj)
	$(C++) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

-include $(dep)   # include all dep files in the makefile

.PHONY: clean
clean:
	rm -f $(obj) $(prefix)eps$(suffix)

.PHONY: cleandep
cleandep:
	rm -f $(dep)
