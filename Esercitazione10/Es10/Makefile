CC =  mpicxx
CFLAGS = -Wall -O3 --std=c++11

a.out: main.o lib.o
	$(CC)  main.o lib.o -o a.out -arch x86_64
main.o: main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS) -arch x86_64
lib.o : lib.cpp lib.h
	$(CC) -c lib.cpp -o lib.o $(CFLAGS) -arch x86_64
clean :
	rm *.o a.out ./dati/*.txt
