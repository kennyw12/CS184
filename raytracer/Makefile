CC = g++
	
RM = /bin/rm -f 
all: main 
main: raytracer.o 
	$(CC) -o raytracer raytracer.o -lpng
raytracer.o: raytracer.cpp
	$(CC) -c raytracer.cpp -o raytracer.o
clean: 
	$(RM) *.o raytracer
 


