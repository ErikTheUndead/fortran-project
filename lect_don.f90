subroutine lect_don

!lol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine pour la lecture des données du problème
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mod_parameter

  IMPLICIT NONE

  select case (Eq)

  case (0) ! Poisson

    open(1,file='initP.data',form='formatted')	

    read(1,*) maillage_type
    read(1,*) Lx,Ly	
    read(1,*) Nx, Ny	
    read(1,*) CL_type
    read(1,*) Neu_D

  case (1) ! Helmholtz

    open(1,file='initH.data',form='formatted')	

    read(1,*) maillage_type	
    read(1,*) Lx,Ly	
    read(1,*) Nx, Ny	
    read(1,*) CL_type
    read(1,*) Neu_D
    read(1,*) dt
    read(1,*) cv
    read(1,*) ech	

  case (2) ! Navier-Stokes Chorin colocatif

    open(1,file='initNS.data',form='formatted')

    read(1,*) maillage_type	
    read(1,*) Lx,Ly	
    read(1,*) Nx, Ny
    read(1,*) cfl
    read(1,*) iter	
    read(1,*) Re
    read(1,*) Nsave
    read(1,*) Temps
    read(1,*) cv

    open(2,file='Nech.data',form='formatted',position='rewind')
    read(2,*) ech
    close(2)

  case (3) ! Navier-Stokes goda colocatif

    open(1,file='initgoda.data',form='formatted')

    read(1,*) maillage_type	
    read(1,*) Lx,Ly	
    read(1,*) Nx, Ny
    read(1,*) cfl
    read(1,*) iter	
    read(1,*) Re
    read(1,*) Nsave
    read(1,*) Temps
    read(1,*) cv

    open(2,file='Nech.data',form='formatted',position='rewind')
    read(2,*) ech
    close(2)

  end select

  close(1)

end subroutine
