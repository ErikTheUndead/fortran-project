Subroutine maillage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine pour la définition du maillage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mod_parameter

  IMPLICIT NONE

  INTEGER :: i,j
  real(kind=8) :: hmin

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Calcul des coordonnées des points du maillage
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!
  !maillage régulier
  !!!!!!!!!!!!!!!!!!

  if ( maillage_type == 0 ) then
  
    do i = 0,Nx 		
      x(i)= real(i)*Lx/real(Nx)
    end do

    do j = 0,Ny
      y(j) = real(j)*Ly/real(Ny)
    end do
    
    Dmin = Lx/Nx

  !!!!!!!!!!!!!!!!!!!!
  !maillage irrégulier
  !!!!!!!!!!!!!!!!!!!!

  else if ( maillage_type == 1) then
  
    do i = 0,Nx 		
      x(i)= 0.5*Lx*(1.-cos(real(i)*Pi/real(Nx)))
    end do

    do j = 0,Ny
      y(j) = 0.5*Ly*(1.-cos(real(j)*Pi/real(Ny)))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calcul de la maille la plus grossière pour le critère CFL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     hmin = 1.
!     
!     do i = 0,Nx-1		  
!       hmin=min(hmin,(x(i+1)-x(i)))		  
!     end do
    
    Dmin = Lx/Nx

  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Calcul des volumes de contrôles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i = 1,Nx-1 ! à l'intérieur
    do j= 1,Ny-1
      Vm(i,j) = (x(i+1)-x(i-1))*(y(j+1)-y(j-1))*0.25
    end do
  end do

  do i = 1,Nx-1 ! Aux bords haut et bas
    Vm(i,0) = (x(i+1)-x(i-1))*(y(1)-y(0))*0.25	
    Vm(i,Ny) = (x(i+1)-x(i-1))*(y(Ny)-y(Ny-1))*0.25
  end do

  do j = 1,Ny-1 ! Aux bords gauche et droit
    Vm(0,j) = (y(j+1)-y(j-1))*(x(1)-x(0))*0.25	
    Vm(Nx,j) = (y(j+1)-y(j-1))*(x(Nx)-x(Nx-1))*0.25
  end do

  ! aux quatre coins

  Vm(0,0) = (y(1)-y(0))*(x(1)-x(0))*0.25	
  Vm(0,Ny) = (y(Ny)-y(Ny-1))*(x(1)-x(0))*0.25
  Vm(Nx,0) = (y(1)-y(0))*(x(Nx)-x(Nx-1))*0.25
  Vm(Nx,Ny) = (y(Ny)-y(Ny-1))*(x(Nx)-x(Nx-1))*0.25

end subroutine maillage	