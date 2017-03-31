BOOST_DIR = /opt/local/include
FFTW_DIR = /opt/local/include
CC	= g++
PARALLEL = -fopenmp #-openmp -openmp-report2
CFLAGS	= -O3 -Wall #-DNDEBUG -pg
LDFLAGS	= -O3 -Wall #-DNDEBUG -pg
INCLUDES = -I$(BOOST_DIR) -I$(HOME)/include
LIBS = -llapack -lblas -L/opt/local/lib -lfftw3 -lboost_system-mt -lboost_filesystem-mt
TARGET	= main 
OBJS	= QuantumDynamics2D.o WaveFunction2D.o main.o

all:	$(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $(PARALLEL) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f $(TARGET) $(OBJS) .nfs* *~ \#* core

.cpp.o:
	$(CC) $(CFLAGS) $(PARALLEL) $(INCLUDES) -c $<
