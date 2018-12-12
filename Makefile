
OBJ = mpiFindMedian.o vp_helper.o main.o
TARGET = vpexec

all: $(OBJ)
	mpicc $^ -o $(TARGET)

%.o: %.c
	mpicc -c $^

clean:
	rm -f *.o *~ $(TARGET)
