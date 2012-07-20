P = main.o emu.o hubble.o
          
LDFLAGS = -O

emu.exe: $(P) $(CSTUFF)
	gcc  -o emu.exe  $(LDFLAGS) $(P) -lgsl -lgslcblas -lm 

