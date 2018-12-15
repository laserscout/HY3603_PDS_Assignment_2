
OBJ = mpiFindMedian.o findMedian.o vp_helper.o main.o
TARGET = vpexec

all: $(OBJ)
	mpicc $^ -o $(TARGET) -lm

%.o: %.c
	mpicc -c $^ -lm

clean:
	rm -f *.o *~ $(TARGET)
