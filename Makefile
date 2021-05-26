sobseq: mod_sobseq.f90 
	gfortran -g $^ -o $@ 

debug: sobseq
	gdb -tui sobseq

test_sobseq: test_sobseq.f90 mod_sobseq.o
	gfortran $^ -o $@

%.o: %.f90
	gfortran -c $<
