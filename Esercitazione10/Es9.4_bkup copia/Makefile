CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o lib.o
	$(CC) main.o lib.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
lib.o : lib.cpp lib.h
	$(CC) -c lib.cpp -o lib.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
