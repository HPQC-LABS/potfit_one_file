ifndef $FC
#  FC = f90             # SUN compilers: f77, f90, f95
   FC = gfortran        # GNU compiler
#  FC = ifort           # Intel compiler   
endif
#                       # Next line: option to choose level of optimization:  make debug=N
ifndef debug
   FFLAGS =  -O2        # debug not specified, so optimize
 else                           
#  FFLAGS =  -q  -C     # debug=1  (basic gtortran debug option)
#  FFLAGS =  -g     -u  # debug=1  (basic debug option)
   ifeq ($(debug),2)
        FFLAGS =  -C -g # debug=2  (higher-level debug option)
    endif
#                       # can add more options for  debug=3 , etc., as desired
  endif
OBJECTS = potfit.o

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o pot.x

#--------------------------------------------------------------
# Sean Mcleod:  13 February 2008
# Nike Dattani:             2008-2025 
#--------------------------------------------------------------
# To run the code, execute:
# ./potfit.x < input.5

# input.5 must have lines 2-3 like the following:
# './2016_04_28.4'             % Name of PotFit input data file
# 'output'                     % Writefile
