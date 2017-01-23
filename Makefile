
all: out.py run.py imclean
	python run.py

imclean:
	rm -f plt/*

stary.dat: dysk.par
	cat $< | diskverta -o stary

out.py: equations.py stary.dat
	python equations.py

clean: imclean
	rm -f *.pyc *.npy
