Geant2txt_float: Geant2txt_float.o 
	gcc -o Geant2txt_float Geant2txt_float.o -lm
Geant2txt_float.o: Geant2txt_float.c
	gcc -c Geant2txt_float.c

clean:
	rm  Geant2txt_float Geant2txt_float.o  -r data -r results -r sino
	mkdir data
	mkdir results
	mkdir sino
# end of Makefile

