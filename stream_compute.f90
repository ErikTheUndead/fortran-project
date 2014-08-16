subroutine stream_compute(u,v)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine pour le calcul de la fonction de courant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mod_parameter
  use mod_thomas
  use mod_inverse

  IMPLICIT NONE

  integer :: i,j 

  REAL(kind=8), allocatable, dimension(:,:) :: u,v

  allocate(u(0:Nx,0:Ny),v(0:Nx,0:Ny))

  ! A l'intérieur
  
  do i=1,Nx-1
    do j=1,Ny-1
      SMs(i,j) = 0.25*(u(i,j+1)-u(i,j-1))*(x(i+1)-x(i-1)) + 0.25*(v(i-1,j)-v(i+1,j))*(y(j+1)-y(j-1))
    end do
  end do	


  ! Terme source à gauche et à droite

  do j=1,Ny-1
    SMs(0,j) = - 0.25*(y(j+1)-y(j-1))*v(1,j)
    SMs(Nx,j) = 0.25*(y(j+1)-y(j-1))*v(Nx-1,j)
  end do

  ! Terme source en bas et en haut

  do i=1,Nx-1
    SMs(i,0) = 0.25*(x(i+1)-x(i-1))*u(i,1)
    SMs(i,Ny) = 0.25*(x(i+1)-x(i-1))*(u(i,Ny)-u(i,Ny-1))
  end do

  ! Sources coins
      
  SMs(0,0) = 0.
  SMs(Nx,0) = 0.

  SMs(0,Ny) = 0. 	  
  SMs(Nx,Ny) = 0. 

  SMs(0,Ny) = 0.25*(x(1)-x(0))*u(0,Ny) 
  SMs(Nx,Ny) = 0.25*(x(Nx)-x(Nx-1))*u(Nx,Ny) 

  call RESOL_THOMAS(Nx,Ny,DgPsi,LwPsi,UpPsi,SMs,Psi)  


! A l'intérieur

!     do i=1,Nx-1
! 	    do j=1,Ny-1
! 	    SMs(i,j) = 0.25*(v(i,j-1)-v(i,j+1))*(x(i+1)-x(i-1)) + 0.25*(u(i+1,j)-u(i-1,j))*(y(j+1)-y(j-1))
! 	    end do
!     end do	  
!   
!   ! Terme source à gauche et à droite
!   
!     do j=1,Ny-1
! 	  SMs(0,j) = 0.25*(y(j+1)-y(j-1))*u(1,j)
! 	  SMs(Nx,j) = - 0.25*(y(j+1)-y(j-1))*u(Nx-1,j)
!     end do
!     
!     ! Terme source en bas et en haut
! 	  
!     do i=1,Nx-1
! 	  SMs(i,0) = - 0.25*(x(i+1)-x(i-1))*v(i,1)
! 	  SMs(i,Ny) = 0.25*(x(i+1)-x(i-1))*v(i,Ny-1)
!     end do
!     
!     ! Sources coins
! 			  
! 	  SMs(0,0) = 0.
! 	  SMs(Nx,0) = 0.
! 	  SMs(0,Ny) = 0. 	  
! 	  SMs(Nx,Ny) = 0.   

end subroutine