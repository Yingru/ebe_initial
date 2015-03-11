
#gfortran -static -ffixed-line-length-none *.f -m32 -O2
#gfortran -static -ffixed-line-length-none *.f -m64 -O2
#gfortran -mcmodel=medium -ffixed-line-length-none *.f -m64 -O2
#gfortran -O3 -ffixed-line-length-none *.f
#gfortran -mcmodel=medium -ffixed-line-length-none -O3 -static -c *.f
gfortran -ffixed-line-length-none -O3 -static -c *.f
gfortran *.o -o initial

echo "initial condition is done"

# done
