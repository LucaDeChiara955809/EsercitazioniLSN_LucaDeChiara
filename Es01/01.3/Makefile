01.3.exe : 01.3.o random.o
	g++ random.o 01.3.o -o 01.3.exe
01.3.o : main.cpp
	g++ -c main.cpp -o 01.3.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o 01.3.exe
