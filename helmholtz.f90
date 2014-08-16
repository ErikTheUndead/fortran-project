subroutine helmholtz

  use mod_parameter
  use mod_thomas
  use mod_inverse

  IMPLICIT NONE

  integer :: i,j,n
  real(kind=8) :: supE,sup2, time1, time2
  real(kind=8), dimension(0:50000) :: supH
  character*50 :: filename
  character*2  :: CL

  allocate(Up2(0:Nx,0:Nx,0:Ny),Dg2(0:Nx,0:Nx,0:Ny),Lw2(0:Nx,0:Nx,0:Ny))
  allocate(f(0:Nx,0:Ny))
  allocate(Tv(0:Nx,0:Ny),T2(0:Nx,0:Ny),T1(0:Nx,0:Ny))
  allocate(Sv(0:Nx,0:Ny),S2(0:Nx,0:Ny))

  cx = (2.*Pi/Lx)
  cy = (2.*Pi/Ly)

  Dg2(:,:,:) = 0.
  Up2(:,:,:) = 0.
  Lw2(:,:,:) = 0.

  supH(:) = 0.
  Tv(:,:) = 0. ! Température sol analytique
  T2(:,:) = 0. ! Température recherchée
  T1(:,:) = 0. ! Température au temps n-1
  Sv(:,:) = 0. ! Source sol analytique
  S2(:,:) = 0. ! Terme source formulé en VF

  ! ====================================================================================
  ! ============ Définition du terme source et température "vraie" =====================
  ! ====================================================================================

  do i = 0,Nx
    do j = 0,Ny
      Sv(i,j) =  - ( cx**2 + cy**2 ) * ( sin(x(i)*cx) * cos(y(j)*cy) )
      Tv(i,j) = sin(x(i)*cx) * cos(y(j)*cy)		
    end do
  end do

  ! ====================================================================================
  ! ====================== Initialisation Dirichlet ====================================
  ! ====================================================================================

  if (CL_type == 0) then

  ! Source en dirichlet = à Température "vraie"

  ! Bords bas et haut

  do i=0,Nx
    Sv(i,0) = Tv(i,0)*Vm(i,0)
    Sv(i,Ny) = Tv(i,Ny)*Vm(i,Ny)	

    Dg(i,i,0) = Vm(i,0)
    Dg(i,i,Ny) = Vm(i,Ny)		
  end do	

  ! Bords gauche et droite

  do j=1,Ny-1
    Sv(0,j) = Tv(0,j)*Vm(0,j)
    Sv(Nx,j) = Tv(Nx,j)*Vm(Nx,j)

    Dg(0,0,j) = Vm(0,j)
    Dg(Nx,Nx,j) = Vm(Nx,j)	
  end do



  ! ====================================================================================
  ! ====================== Initialisation Neumann ======================================
  ! ====================================================================================

  else if (CL_type == 1) then
  
  !!!!!!!!!!!!!!!!!!!
  ! Définition de f : 
  !!!!!!!!!!!!!!!!!!!
  
  ! parois haute et basse

  do i=1,Nx-1
    f(i,0) = cy*sin(x(i)*cx)*sin(y(0)*cy)
    f(i,Ny) = - cy*sin(x(i)*cx)*sin(y(Ny)*cy)
  end do

  ! parois gauche et droite

  do j=1,Ny-1
    f(0,j) = - cx*cos(x(0)*cx)*cos(y(j)*cy)
    f(Nx,j) = cx*cos(x(Nx)*cx)*cos(y(j)*cy)
  end do

  ! Bords gauche et droit

  do j=1,Ny-1

    Sv(0,j) = Sv(0,j)*Vm(0,j) - f(0,j)*(y(j+1)-y(j-1))*0.5		
    Sv(Nx,j) = Sv(Nx,j)*Vm(Nx,j) - f(Nx,j)*(y(j+1)-y(j-1))*0.5

    Dg(0,0,j) = - 0.5*(x(1)-x(0))/(y(j)-y(j-1)) - 0.5*(y(j+1)-y(j-1))/(x(1)-x(0)) - 0.5*(x(1)-x(0))/(y(j+1)-y(j)) ! T(0,j)
    Dg(0,1,j) = 0.5*(y(j+1)-y(j-1))/(x(1)-x(0)) ! T(1,j)
    Lw(0,0,j) = 0.5*(x(1)-x(0))/(y(j)-y(j-1)) ! T(0,j-1)
    Up(0,0,j) = 0.5*(x(1)-x(0))/(y(j+1)-y(j)) ! T(0,j+1)

    Dg(Nx,Nx,j) = - 0.5*(x(Nx)-x(Nx-1))/(y(j)-y(j-1)) - 0.5*(y(j+1)-y(j-1))/(x(Nx)-x(Nx-1)) - 0.5*(x(Nx)-x(Nx-1))/(y(j+1)-y(j)) ! T(Nx,j)
    Dg(Nx,Nx-1,j) = 0.5*(y(j+1)-y(j-1))/(x(Nx)-x(Nx-1)) ! T(Nx-1,j)
    Lw(Nx,Nx,j) = 0.5*(x(Nx)-x(Nx-1))/(y(j)-y(j-1)) ! T(Nx,j-1)
    Up(Nx,Nx,j) = 0.5*(x(Nx)-x(Nx-1))/(y(j+1)-y(j)) ! T(Nx,j+1)

  end do		

  ! Bords haut et bas

  do i=1,Nx-1

    Sv(i,Ny) = Sv(i,Ny)*Vm(i,Ny) - f(i,Ny)*(x(i+1)-x(i-1))*0.5		
    Sv(i,0) = Sv(i,0)*Vm(i,0) - f(i,0)*(x(i+1)-x(i-1))*0.5

    Dg(i,i,0) = - 0.5*(y(1)-y(0))/(x(i)-x(i-1)) - 0.5*(x(i+1)-x(i-1))/(y(1)-y(0)) - 0.5*(y(1)-y(0))/(x(i+1)-x(i)) ! T(i,0)
    Dg(i,i-1,0) = 0.5*(y(1)-y(0))/(x(i)-x(i-1)) ! T(i-1,0)
    Dg(i,i+1,0) = 0.5*(y(1)-y(0))/(x(i+1)-x(i)) ! T(i+1,0)
    Up(i,i,0) = 0.5*(x(i+1)-x(i-1))/(y(1)-y(0)) ! T(i,1)

    Dg(i,i,Ny) = - 0.5*(y(Ny)-y(Ny-1))/(x(i)-x(i-1)) - 0.5*(x(i+1)-x(i-1))/(y(Ny)-y(Ny-1)) - 0.5*(y(Ny)-y(Ny-1))/(x(i+1)-x(i)) ! T(i,Ny)
    Dg(i,i-1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(i)-x(i-1)) ! T(i-1,0)
    Dg(i,i+1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(i+1)-x(i)) ! T(i+1,0)
    Lw(i,i,Ny) = 0.5*(x(i+1)-x(i-1))/(y(Ny)-y(Ny-1)) ! T(i,Ny-1)

  end do


  if (Neu_D == 1) then ! 4 points dirichlets aux 4 coins

    Sv(Nx,Ny) = Tv(Nx,Ny)*Vm(Nx,Ny)
    Dg(Nx,Nx,Ny) = Vm(Nx,Ny)

    Sv(0,Ny) = Tv(0,Ny)*Vm(0,Ny)
    Dg(0,0,Ny) = Vm(0,Ny)

    Sv(Nx,0) = Tv(Nx,0)*Vm(Nx,0)
    Dg(Nx,Nx,0) = Vm(Nx,0)

    Sv(0,0) = Tv(0,0)*Vm(0,0)
    Dg(0,0,0) = Vm(0,0)

  else  ! 1 point dirichlet au coin Nx,Ny

  ! Point (0,0)		

    Sv(0,0) = Sv(0,0)*Vm(0,0) - f(1,0)*(x(1)-x(0))*0.5 - f(0,1)*(y(1)-y(0))*0.5

    Dg(0,0,0) = - 0.5*(y(1)-y(0))/(x(1)-x(0)) - 0.5*(x(1)-x(0))/(y(1)-y(0)) ! T(0,0)
    Dg(0,1,0) = 0.5*(y(1)-y(0))/(x(1)-x(0)) ! T(1,0)
    Up(0,0,0) = 0.5*(x(1)-x(0))/(y(1)-y(0)) ! T(0,1)

    ! Point (0,Ny) 		

    Sv(0,Ny) = Sv(0,Ny)*Vm(0,Ny) - f(1,Ny)*(x(1)-x(0))*0.5 - f(0,Ny-1)*(y(Ny)-y(Ny-1))*0.5

    Dg(0,0,Ny) = - 0.5*(y(Ny)-y(Ny-1))/(x(1)-x(0)) - 0.5*(x(1)-x(0))/(y(Ny)-y(Ny-1)) ! T(0,Ny)
    Dg(0,1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(1)-x(0)) ! T(1,Ny)
    Lw(0,0,Ny) = 0.5*(x(1)-x(0))/(y(Ny)-y(Ny-1)) ! T(0,Ny-1)

    ! Point (Nx,0) 		

    Sv(Nx,0) = Sv(Nx,0)*Vm(Nx,0) - f(Nx-1,0)*(x(Nx)-x(Nx-1))*0.5 - f(Nx,1)*(y(1)-y(0))*0.5

    Dg(Nx,Nx,0) = - 0.5*(y(1)-y(0))/(x(Nx)-x(Nx-1)) - 0.5*(x(Nx)-x(Nx-1))/(y(1)-y(0)) ! T(Nx,0)
    Dg(Nx,Nx-1,0) = 0.5*(y(1)-y(0))/(x(Nx)-x(Nx-1)) ! T(Nx-1,0) 
    Up(Nx,Nx,0) = 0.5*(x(Nx)-x(Nx-1))/(y(1)-y(0)) ! T(Nx,1)

    ! Point (Nx,Ny) Point Dirichlet

    Sv(Nx,Ny) = Tv(Nx,Ny)*Vm(Nx,Ny)
    Dg(Nx,Nx,Ny) = Vm(Nx,Ny)

  end if

  end if	

  ! ====================================================================================
  ! ================= Réecriture du terme source en VF =================================
  ! ====================================================================================

  ! A l'intérieur 

  do i = 1,Nx-1
    do j = 1,Ny-1
      Sv(i,j) = Sv(i,j)*Vm(i,j)
    end do
  end do

  S2(:,:) = Sv(:,:)

  ! on souhaite conserver Up,Dg et Lw à part.

  Up2(:,:,:) = Up(:,:,:)
  Lw2(:,:,:) = Lw(:,:,:)
  Dg2(:,:,:) = Dg(:,:,:)

  ! ====================================================================================
  ! ============ Ecriture du problème de Helmholtz et résolution =======================
  ! ====================================================================================

  ! ========================= réecriture des matrices ==================================

  do j = 1,Ny-1

  ! à l'intérieur

  do i = 1,Nx-1

    Dg2(i,i,j) = Dg(i,i,j)*dt - Vm(i,j)

    Dg2(i,i+1,j) = Dg(i,i+1,j)*dt
    Dg2(i,i-1,j) = Dg(i,i-1,j)*dt

    Lw2(i,i,j) = Lw(i,i,j)*dt
    Up2(i,i,j) = Up(i,i,j)*dt

  end do

  end do

  ! ===========================================================================

  call INIT_THOMAS(Nx,Ny,Dg2,Lw2,Up2)

  n = 1
  time1 = 0.
  time2 = 0.	

  ! ============= début de la boucle pour le problème de Helmholtz ============	


  call cpu_time(time1)

  do while ( n < 10 .or. abs(supH(n)-supH(n-1)) > cv )

    do i = 1,Nx-1
      do j = 1,Ny-1
	S2(i,j) = Sv(i,j)*dt - T1(i,j)*Vm(i,j)			
      end do
    end do

    call RESOL_THOMAS(Nx,Ny,Dg2,Lw2,Up2,S2,T2)

    do i=0,Nx
      do j=0,Ny
	supH(n) = max(supH(n),abs(T2(i,j)-T1(i,j)))
      enddo
    enddo

    T1(:,:) = T2(:,:)
    n = n+1

  end do

  ! ============= Ecriture dans un fichier =====================================

  supE = 0.

  ! Erreur max à convergence

  do i=0,Nx
    do j=0,Ny
      supE = max(supE,abs(T2(i,j)-Tv(i,j)))
    enddo
  enddo

  call cpu_time(time2)

  if (CL_type == 0) then
    write(CL,'(a)') 'Di'
  else if (CL_type == 1) then
    write(CL,'(a)') 'Ne'
  end if 

  ! finesse de la plus grande maille

  do i=1,Nx-1				
    sup2= max(sup2,abs(x(i)-x(i-1)))	
  end do	

  write(filename,'(a,a,a,i1.1,a,i2.2,a)')'ConvergenceHelmholtz_',CL,'_',maillage_type,'_echantillon_',ech,'_dt.dat'
  open(1,file=filename,position="append",form='formatted',status='unknown')
  write(1,*) dt,supE
  close(1)

  write(*,*) n,sup2,Nx,Ny,time2-time1,supE

  deallocate(f)
  deallocate(Tv,T2,T1)
  deallocate(Sv,S2)

end subroutine

