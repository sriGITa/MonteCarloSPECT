# This is Makefile to compile Jaszczak Phantom 

Jaszczak: Jaszczak.o 
	gcc -o Jaszczak Jaszczak.o -lm

Jaszczak.o: Jaszczak.c
	gcc -c Jaszczak.c
# improterty depends on one source files

clean:
	rm  Jaszczak Jaszczak.o -r data
	mkdir data
# end of Makefile

