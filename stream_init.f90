subroutine stream_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine pour initialiser les matrices D,L et U pour le calcul de la fonction de courant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mod_parameter
  use mod_thomas
  use mod_inverse

  IMPLICIT NONE

  integer :: i,j      

  DgPsi = Dg
  UpPsi = Up
  LwPsi = Lw

  !!!!!!!!!!!!!!!!!!!!!!!
  !Bords gauche et droite
  !!!!!!!!!!!!!!!!!!!!!!!

  do j=1,Ny-1		

    DgPsi(0,1,j) = 0.5*(y(j+1)-y(j-1))/(x(1)-x(0)) ! Psi(1,j)
    LwPsi(0,0,j) = 0.5*(x(1)-x(0))/(y(j)-y(j-1)) ! Psi(0,j-1)
    UpPsi(0,0,j) = 0.5*(x(1)-x(0))/(y(j+1)-y(j)) ! Psi(0,j+1)
    DgPsi(0,0,j) = - DgPsi(0,1,j) - LwPsi(0,0,j) - UpPsi(0,0,j) ! Psi(0,j)

    DgPsi(Nx,Nx-1,j) = 0.5*(y(j+1)-y(j-1))/(x(Nx)-x(Nx-1)) ! Psi(Nx-1,j)
    LwPsi(Nx,Nx,j) = 0.5*(x(Nx)-x(Nx-1))/(y(j)-y(j-1)) ! Psi(Nx,j-1)
    UpPsi(Nx,Nx,j) = 0.5*(x(Nx)-x(Nx-1))/(y(j+1)-y(j)) ! Psi(Nx,j+1)
    DgPsi(Nx,Nx,j) = - DgPsi(Nx,Nx-1,j) - LwPsi(Nx,Nx,j) - UpPsi(Nx,Nx,j) ! Psi(Nx,j)

  end do	


  !!!!!!!!!!!!!!!!!!
  !Bords haut et bas
  !!!!!!!!!!!!!!!!!!

  do i=1,Nx-1

    DgPsi(i,i-1,0) = 0.5*(y(1)-y(0))/(x(i)-x(i-1)) ! Psi(i-1,0)
    DgPsi(i,i+1,0) = 0.5*(y(1)-y(0))/(x(i+1)-x(i)) ! Psi(i+1,0)
    UpPsi(i,i,0) = 0.5*(x(i+1)-x(i-1))/(y(1)-y(0)) ! Psi(i,1)
    DgPsi(i,i,0) = - DgPsi(i,i-1,0) - DgPsi(i,i+1,0) - UpPsi(i,i,0)! Psi(i,0)

    DgPsi(i,i-1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(i)-x(i-1)) ! Psi(i-1,Ny)
    DgPsi(i,i+1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(i+1)-x(i)) ! Psi(i+1,Ny)
    LwPsi(i,i,Ny) = 0.5*(x(i+1)-x(i-1))/(y(Ny)-y(Ny-1)) ! Psi(i,Ny-1)
    DgPsi(i,i,Ny) = - DgPsi(i,i-1,Ny) - DgPsi(i,i+1,Ny) - LwPsi(i,i,Ny)  ! Psi(i,Ny)

  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Conditions limites aux coins : 1 coin Dirichlet + 3 coins Neumann	
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Point (Nx,0) Neumann
	  
  DgPsi(Nx,Nx-1,0) = 0.5*(y(1)-y(0))/(x(Nx)-x(Nx-1)) 
  UpPsi(Nx,Nx,0) = 0.5*(x(Nx)-x(Nx-1))/(y(1)-y(0))
  DgPsi(Nx,Nx,0) = - DgPsi(Nx,Nx-1,0) - UpPsi(Nx,Nx,0)

  !Point (0,0) Neumann

  DgPsi(0,1,0) = 0.5*(y(1)-y(0))/(x(1)-x(0)) 
  UpPsi(0,0,0) = 0.5*(x(1)-x(0))/(y(1)-y(0))
  DgPsi(0,0,0) = - DgPsi(0,1,0) - UpPsi(0,0,0)

  !Point (0,Ny) Dirichlet

  DgPsi(0,0,Ny) = 1.

  !Point (Nx,Ny) Neumann 	

  DgPsi(Nx,Nx-1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(Nx)-x(Nx-1)) 
  LwPsi(Nx,Nx,Ny) = 0.5*(x(Nx)-x(Nx-1))/(y(Ny)-y(Ny-1))
  DgPsi(Nx,Nx,Ny) = - DgPsi(Nx,Nx-1,Ny) - LwPsi(Nx,Nx,Ny)

  call INIT_THOMAS(Nx,Ny,DgPsi,LwPsi,UpPsi)


end subroutine
