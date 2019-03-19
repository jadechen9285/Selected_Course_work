OBJECTS = LinAl.o Driver_HH.o

MODULES = LinAl.mod

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe

%.o: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

%.mod: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl.exe
