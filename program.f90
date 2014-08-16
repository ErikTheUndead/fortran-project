program Program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Programme principale pour la résolution des différents problèmes de Poisson, Helmholtz et Navier-Stokes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mod_parameter !module regroupant toute les variables communes aux différentes subroutines
  use mod_thomas
  use mod_inverse
  use mod_save_vtk 

  IMPLICIT NONE

  INTEGER :: i,j

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Affichage d'un menu dans la console pour le choix de la résolution
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1 write(*,*) 'Choisissez le type de problème à résoudre :'
  write(*,*) 'Résolution conduction : Tapez 0 pour Poisson, 1 pour Helmholtz'
  write(*,*) 'Résolution NS : Tapez 2 pour Chorin colocatif, 3 pour goda colocatif'

  read(*,*) Eq	

  select case (Eq)
    case(0)
      write(*,*) 'Choix résolution Poisson (initP.data)'
    case(1)
      write(*,*) 'Choix résolution Helmholtz (initH.data)'
    case(2)
      write(*,*) 'Choix résolution Chorin colocatif (initNS.data)'
    case(3)
      write(*,*) 'Choix résolution Goda colocatif (initgoda.data)'
    case default
      write(*,*) 'Choix invalide, retour au menu'	
      go to 1
  end select

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
  !Lecture des paramètres "lect_don.f90" dans des fichiers, selon le type de résolution voulu.
  !initNS.data (Chorin), initgoda.data (goda), initP.data (Poisson) & initH.data (Helmholtz)!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call lect_don

  allocate(x(0:Nx))
  allocate(y(0:Ny))
  allocate(Vm(0:Nx,0:Ny))	

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialisation du maillage "maillage.f90"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call maillage

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
  !Ecriture du laplacien de réference sans les conditions limites
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(Up(0:Nx,0:Nx,0:Ny),Dg(0:Nx,0:Nx,0:Ny),Lw(0:Nx,0:Nx,0:Ny))

  Dg(:,:,:) = 0.
  Up(:,:,:) = 0.
  Lw(:,:,:) = 0.

  do j=1,Ny-1
    do i=1,Nx-1
      Dg(i,i-1,j) = 0.5*(y(j+1)-y(j-1))/(x(i)-x(i-1))
      Dg(i,i+1,j) = 0.5*(y(j+1)-y(j-1))/(x(i+1)-x(i))
      Lw(i,i,j) =  0.5*(x(i+1)-x(i-1))/(y(j)-y(j-1))
      Up(i,i,j) = 0.5*(x(i+1)-x(i-1))/(y(j+1)-y(j))
      Dg(i,i,j) = -(Dg(i,i-1,j)+Dg(i,i+1,j)+Lw(i,i,j)+Up(i,i,j))
    end do	
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
  ! Choix des subroutines pour la résolution
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  select case (Eq)
    case(0)
      call poisson
    case(1)
      call helmholtz
    case(2)
      call init
      call NS
    case(3)
      call init
      call goda
  end select	


  deallocate(x)
  deallocate(y)
  deallocate(Vm)
  deallocate(Up,Dg,Lw)


END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine pour l'initialisation de l'algorithme de Chorin et Goda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init

  use mod_parameter

  implicit none  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !la signification des variables est décrite dans mod_parameter
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(usx(0:Nx,0:Ny),usy(0:Nx,0:Ny))
  allocate(ux1(0:Nx,0:Ny),uy1(0:Nx,0:Ny))
  allocate(ux2(0:Nx,0:Ny),uy2(0:Nx,0:Ny))
  allocate(p(0:Nx,0:Ny),phi(0:Nx,0:Ny),p1(0:Nx,0:Ny),p2(0:Nx,0:Ny))
  allocate(SMx(0:Nx,0:Ny),SMy(0:Nx,0:Ny),SMp(0:Nx,0:Ny))

  allocate(Up2(0:Nx,0:Nx,0:Ny),Dg2(0:Nx,0:Nx,0:Ny),Lw2(0:Nx,0:Nx,0:Ny))
  allocate(Up3(0:Nx,0:Nx,0:Ny),Dg3(0:Nx,0:Nx,0:Ny),Lw3(0:Nx,0:Nx,0:Ny))  
  allocate(DgPsi(0:Nx,0:Nx,0:Ny),UpPsi(0:Nx,0:Nx,0:Ny),LwPsi(0:Nx,0:Nx,0:Ny))

  allocate(Psi(0:Nx,0:Ny))
  allocate(SMs(0:Nx,0:Ny))

  allocate(sup(0:iter,1:3))

  !U0 = Re*10.0**(-3)/Lx
  U0 = 1.

  write(*,*)
  write(*,*) '======================================================'
  write(*,*) 'U0 = ', U0
  write(*,*) 'Re = ', Re
  dt = cfl*Dmin/abs(U0)
  write(*,*) 'dt =', dt,'(s)'
  write(*,*) 'Temps simulation max : ', Temps,' s'
  write(*,*) 'convergence :', cv
  write(*,*) '======================================================'
  write(*,*)

end subroutine
