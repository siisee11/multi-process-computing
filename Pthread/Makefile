CC = g++
CFLAGS = -W -Wall
TARGET = multithreadcpp
OBJECTS = project1.o

all : $(TARGET)

$(TARGET) : $(OBJECTS) 
	$(CC) $(CFLAGS) -o $@ $^ -lpthread -O2

clean :
	rm *.o multithread
