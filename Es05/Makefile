05.1.exe : 05.1.o random.o
	g++ Posizione.o random.o 05.1.o -o 05.1.exe
05.1.o : main.cpp
	g++ -c main.cpp -o 05.1.o
Posizione.o : Posizione.cpp Posizione.h
	g++ -c Posizione.cpp -o Posizione.o 
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o 05.1.exe
