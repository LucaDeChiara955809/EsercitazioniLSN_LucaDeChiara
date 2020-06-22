10.2.exe : 10.2.o random.o
	g++ random.o Individuo.o 10.2.o -o 10.2.exe
10.2.o : main.cpp
	g++ -c main.cpp -o 10.2.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
Individuo.o : Individuo.cpp Individuo.h
	g++ -c Individuo.cpp -o Individuo.o 
clean :
	rm *.o 10.2.exe
