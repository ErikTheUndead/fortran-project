FC = gfortran
LD = gfortran
# 
FFLAGS = -C -fbounds-check 
# 
LDFLAGS = -fbounds-check 
#  
LINK = -c


.SUFFIXES : .f90 .o

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90

OBJET = mod_parameter.o maillage.o lect_don.o mod_inverse.o mod_thomas.o mod_save_vtk.o poisson.o helmholtz.o NS.o NS_s.o goda.o stream_init.o stream_compute.o


program.x : program.o $(OBJET)
	$(LD) $(LDFLAGS) -o program.x program.o $(OBJET) $(LIBS)

mod_parameter.o : mod_parameter.f90
	$(FC) $(FFLAGS) $(LINK) mod_parameter.f90

mod_thomas.o : mod_thomas.f90 
	$(FC) $(FFLAGS) $(LINK) mod_thomas.f90

mod_inverse.o : mod_inverse.f90 
	$(FC) $(FFLAGS) $(LINK) mod_inverse.f90

mod_save_vtk.o : mod_save_vtk.f90 
	$(FC) $(FFLAGS) $(LINK) mod_save_vtk.f90
	

program.o : program.f90 mod_parameter.o lect_don.o mod_inverse.o mod_thomas.o mod_save_vtk.o poisson.o helmholtz.o NS.o NS_s.o goda.o stream_init.o stream_compute.o
	$(FC) $(FFLAGS) $(LINK) program.f90

lect_don.o: lect_don.f90 mod_parameter.o
	$(FC) $(FFLAGS) $(LINK) lect_don.f90

maillage.o: maillage.f90 mod_parameter.o
	$(FC) $(FFLAGS) $(LINK) maillage.f90

poisson.o: poisson.f90 mod_parameter.o  mod_inverse.o mod_thomas.o
	$(FC) $(FFLAGS) $(LINK) poisson.f90

helmholtz.o: helmholtz.f90 mod_parameter.o mod_inverse.o mod_thomas.o
	$(FC) $(FFLAGS) $(LINK) helmholtz.f90  

NS.o: NS.f90 mod_parameter.o mod_inverse.o mod_thomas.o mod_save_vtk.o 
	$(FC) $(FFLAGS) $(LINK) NS.f90

NS_s.o: NS_s.f90 mod_parameter.o mod_inverse.o mod_thomas.o mod_save_vtk.o
	$(FC) $(FFLAGS) $(LINK) NS_s.f90

goda.o: goda.f90 mod_parameter.o mod_inverse.o mod_thomas.o mod_save_vtk.o 
	$(FC) $(FFLAGS) $(LINK) goda.f90
	
stream_init.o: stream_init.f90 mod_parameter.o mod_inverse.o mod_thomas.o 
	$(FC) $(FFLAGS) $(LINK) stream_init.f90
	
stream_compute.o: stream_compute.f90 mod_parameter.o mod_inverse.o mod_thomas.o 
	$(FC) $(FFLAGS) $(LINK) stream_compute.f90
	
clean:
	rm -f *.o *.a *.x *.mod

cleand: 
	rm -f *.txt *.vtr *.data~ *.f90~ *.dat~

again: clean all




