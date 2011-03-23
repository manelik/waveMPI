



MPIFCOMP = mpif90

axionda:	axionda.f90
		$(MPIFCOMP) axionda.f90 -o axionda

onda:		onda.f90
		$(MPIFCOMP) onda.f90 -o onda


clean:		
		rm *.rl *~ 
 
