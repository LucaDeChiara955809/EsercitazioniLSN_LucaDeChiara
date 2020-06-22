01.2.exe : 01.2.o random.o
	g++ random.o 01.2.o -o 01.2.exe
01.2.o : main.cpp
	g++ -c main.cpp -o 01.2.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o 01.2.exe
