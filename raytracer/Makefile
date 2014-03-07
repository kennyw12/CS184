CC = g++
CFLAGS = -g
LDFLAGS = -lglut -lGLU

	
RM = /bin/rm -f 
all: main 
main: raytracer.o 
	$(CC) $(CFLAGS) -o as2 raytracer.o $(LDFLAGS) 
raytracer.o: raytracer.cpp
	$(CC) $(CFLAGS) -c raytracer.cpp -o raytracer.o
clean: 
	$(RM) *.o as2