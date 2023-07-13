CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.x : main.o random.o
	$(CC) random.o main.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random/random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.x
clean_H :
	rm ./H/pos.out ./H/result.out
clean_SA :
	rm ./SA/*.out ./SA/*.txt
	
