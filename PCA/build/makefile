# Linux
CC=g++
LIBS= -std=c++11
OBJS= main.o

# Compiler flags
CXXFLAGS= -g -Wextra

all: ../bin/Param

../bin/Param: $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o ../bin/Param $(LIBS)

main.o: ../src/main.cpp ../src/color.h ../src/linear.h
	$(CC) $(CXXFLAGS) -c ../src/main.cpp $(LIBS)

clean:
	rm -rf ../bin/Param ./*~ *~ *.o
