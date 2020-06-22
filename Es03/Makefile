03.1.exe : 03.1.o random.o
	g++ random.o 03.1.o -o 03.1.exe
03.1.o : main.cpp
	g++ -c main.cpp -o 03.1.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o 03.1.exe
