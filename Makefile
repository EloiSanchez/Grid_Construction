comp = gfortran
opt = -Wall
#opt = -O3

grid_constr :
	$(comp) $(opt) grid_constr.f90 -o grid_constr.x

.PHONY: clean
clean :
	rm -rf *.x
