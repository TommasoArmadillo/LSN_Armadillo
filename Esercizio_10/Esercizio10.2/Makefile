CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.cpp main.h random.o city.o cromosome.o
	mpic++ random.o main.cpp city.o cromosome.o -o main.exe
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)
cromosome.o : cromosome.cpp cromosome.h
	$(CC) -c cromosome.cpp -o cromosome.o $(CFLAGS)
clean :
	rm *.o main.exe
