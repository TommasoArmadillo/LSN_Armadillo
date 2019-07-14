CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o city.o cromosome.o
	$(CC) random.o main.o city.o cromosome.o -o main.exe
main.o : main.cpp main.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)
cromosome.o : cromosome.cpp cromosome.h
	$(CC) -c cromosome.cpp -o cromosome.o $(CFLAGS)
clean :
	rm *.o main.exe
