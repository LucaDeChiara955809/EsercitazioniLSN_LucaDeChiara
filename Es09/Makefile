09.1.exe : 09.1.o random.o
	g++ random.o Individuo.o 09.1.o -o 09.1.exe
09.1.o : main.cpp
	g++ -c main.cpp -o 09.1.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
Individuo.o : Individuo.cpp Individuo.h
	g++ -c Individuo.cpp -o Individuo.o 
clean :
	rm *.o 09.1.exe
