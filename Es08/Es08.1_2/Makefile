08.1_2.exe : 08.1_2.o random.o
	g++ random.o 01.1.o -o 01.1.exe
08.1_2.o : main.cpp
	g++ -c main.cpp -o 01.1.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o 08.1_2.exe
