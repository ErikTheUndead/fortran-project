subroutine NS

  use mod_parameter
  use mod_thomas
  use mod_inverse
  use mod_save_vtk

  IMPLICIT NONE

  integer :: i,j,k,n,inc,inc2
  real(kind=8) :: imax,jmax
  real(kind=8) :: sup2,sup3
  real(kind=8) :: tps1, tps2 

  Dg2 = 0.
  Up2 = 0.
  Lw2 = 0.
  Dg3 = 0.
  Up3 = 0.
  Lw3 = 0.	

  Dg2 = Dg
  Up2 = Up
  Lw2 = Lw
  
  Psi = 0.

  ! ====================================================================================
  ! ============ Ecriture du problème de Navier-Stokes =================================
  ! ====================================================================================

  ! ========================= Initialisation ===================================
  ! ========================= réecriture des matrices ==========================
  ! ========================= pour le calcul de u* et v* =======================

  do j = 1,Ny-1
    do i = 1,Nx-1
      Dg2(i,i,j) = Vm(i,j) - Dg2(i,i,j)*dt/Re 
      Dg2(i,i+1,j) = -Dg2(i,i+1,j)*dt/Re
      Dg2(i,i-1,j) = -Dg2(i,i-1,j)*dt/Re
      Lw2(i,i,j) = -Lw2(i,i,j)*dt/Re
      Up2(i,i,j) = -Up2(i,i,j)*dt/Re
    end do
  end do

  !========================== aux CLs =========================================

  ! CL dirichlet pour la vitesse

  do i=0,Nx
    Dg2(i,i,0) = Vm(i,0)
    Dg2(i,i,Ny) = Vm(i,Ny)
    Lw2(i,i,0) = 0.
    Lw2(i,i,Ny) = 0.	
    Up2(i,i,0) = 0.
    Up2(i,i,Ny) = 0.
  end do	

  do j=1,Ny-1
    Dg2(0,0,j) = Vm(0,j)
    Dg2(Nx,Nx,j) = Vm(Nx,j)
    Lw2(0,0,j) = 0.
    Lw2(Nx,Nx,j) = 0.
    Up2(0,0,j) = 0.
    Up2(Nx,Nx,j) = 0.
  end do	

  ! =================== CLs Dirichlet Second membre ============================= 

  ! Bords haut et bas

  do i=0,Nx
    SMx(i,0) = 0.
    SMx(i,Ny) = U0*Vm(i,Ny)	
  end do

  ! Bords gauche et droite

  do j=1,Ny-1
    SMx(0,j) = 0.
    SMx(Nx,j) = 0.
  end do

  ! projeté en y : Nul partout

  ! Bords haut et bas

  do i=0,Nx
    SMy(i,0) = 0.
    SMy(i,Ny) = 0.	
  end do

  ! Bords gauche et droite

  do j=1,Ny-1
    SMy(0,j) = 0.
    SMy(Nx,j) = 0.
  end do

  call INIT_THOMAS(Nx,Ny,Dg2,Lw2,Up2)

  ! ===================== Conditions initiales ================================

  uy1 = 0.
  ux1 = 0.
  sup = 0.

  do i=0,Nx
    ux1(i,Ny) = U0
  end do	

  ! ========================= réecriture des matrices ==========================
  ! ========================= pour le calcul P(n+1) ============================

  Up3 = Up
  Lw3 = Lw
  Dg3 = Dg

  ! ========================= CL Neumann =======================================

  ! Bords gauche et droite

  do j=1,Ny-1		

    Dg3(0,1,j) = 0.5*(y(j+1)-y(j-1))/(x(1)-x(0)) ! p(1,j)
    Lw3(0,0,j) = 0.5*(x(1)-x(0))/(y(j)-y(j-1)) ! p(0,j-1)
    Up3(0,0,j) = 0.5*(x(1)-x(0))/(y(j+1)-y(j)) ! p(0,j+1)
    Dg3(0,0,j) = - Dg3(0,1,j) - Lw3(0,0,j) - Up3(0,0,j) ! p(0,j)

    Dg3(Nx,Nx-1,j) = 0.5*(y(j+1)-y(j-1))/(x(Nx)-x(Nx-1)) ! p(Nx-1,j)
    Lw3(Nx,Nx,j) = 0.5*(x(Nx)-x(Nx-1))/(y(j)-y(j-1)) ! p(Nx,j-1)
    Up3(Nx,Nx,j) = 0.5*(x(Nx)-x(Nx-1))/(y(j+1)-y(j)) ! p(Nx,j+1)
    Dg3(Nx,Nx,j) = - Dg3(Nx,Nx-1,j) - Lw3(Nx,Nx,j) - Up3(Nx,Nx,j) ! p(Nx,j)

  end do		

  ! Bords haut et bas

  do i=1,Nx-1

    Dg3(i,i-1,0) = 0.5*(y(1)-y(0))/(x(i)-x(i-1)) ! p(i-1,0)
    Dg3(i,i+1,0) = 0.5*(y(1)-y(0))/(x(i+1)-x(i)) ! p(i+1,0)
    Up3(i,i,0) = 0.5*(x(i+1)-x(i-1))/(y(1)-y(0)) ! p(i,1)
    Dg3(i,i,0) = - Dg3(i,i-1,0) - Dg3(i,i+1,0) - Up3(i,i,0)! p(i,0)

    Dg3(i,i-1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(i)-x(i-1)) ! p(i-1,Ny)
    Dg3(i,i+1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(i+1)-x(i)) ! p(i+1,Ny)
    Lw3(i,i,Ny) = 0.5*(x(i+1)-x(i-1))/(y(Ny)-y(Ny-1)) ! p(i,Ny-1)
    Dg3(i,i,Ny) = - Dg3(i,i-1,Ny) - Dg3(i,i+1,Ny) - Lw3(i,i,Ny)  ! p(i,Ny)

  end do

  ! 1 coin Dirichlet pour la pression + 3 coins Neumann


  !Point (Nx,0) Dirichlet

  Dg3(Nx,Nx,0) = 1.
  SMp(Nx,0) = 0.	

  !Point (0,0) Neumann

  Dg3(0,1,0) = 0.5*(y(1)-y(0))/(x(1)-x(0)) 
  Up3(0,0,0) = 0.5*(x(1)-x(0))/(y(1)-y(0))
  Dg3(0,0,0) = - Dg3(0,1,0) - Up3(0,0,0)

  ! Point (0,Ny) Neumann	

  Dg3(0,1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(1)-x(0)) 
  Lw3(0,0,Ny) = 0.5*(x(1)-x(0))/(y(Ny)-y(Ny-1)) 
  Dg3(0,0,Ny) = - Dg3(0,1,Ny) - Lw3(0,0,Ny)

  ! Point (Nx,Ny) Neumann 	

  Dg3(Nx,Nx-1,Ny) = 0.5*(y(Ny)-y(Ny-1))/(x(Nx)-x(Nx-1)) 
  Lw3(Nx,Nx,Ny) = 0.5*(x(Nx)-x(Nx-1))/(y(Ny)-y(Ny-1))
  Dg3(Nx,Nx,Ny) = - Dg3(Nx,Nx-1,Ny) - Lw3(Nx,Nx,Ny)

  ! CL pour coin Neumann

  SMp(0,0) = 0.
  SMp(0,Ny) = 0. 
  SMp(Nx,Ny) = 0.

  call INIT_THOMAS(Nx,Ny,Dg3,Lw3,Up3)

  inc = 0
  inc2 = 0
  
  ! Initialisation des matrices pour la fonction de courant

  call stream_init

  ! ============================================================================
  ! ========================= Calcul ===========================================
  ! ============================================================================

  n = 1

  sup(0,1) = 0.
  sup(1,1) = 1.
  
  call cpu_time(tps1)

  do while ( abs(sup(n,1)) > cv .and. (n-1)*dt < Temps)

  n=n+1	

    usx = 0.
    usy = 0.

    ! ========================= Etape 1 : prédiction de u* ======================
    ! ===========================================================================

    ! =================== Ecriture du Second Membre =============================
    ! =================== à l'intérieur =========================================

    do i=1,Nx-1
      do j=1,Ny-1		
	SMx(i,j) = - (0.25*(y(j+1)-y(j-1))*(ux1(i+1,j)**2-ux1(i-1,j)**2) + 0.25*(ux1(i,j+1)*uy1(i,j+1) - &
	& ux1(i,j-1)*uy1(i,j-1))*(x(i+1)-x(i-1)))*dt + ux1(i,j)*Vm(i,j)
	SMy(i,j) = - (0.25*(y(j+1)-y(j-1))*(ux1(i+1,j)*uy1(i+1,j) - ux1(i-1,j)*uy1(i-1,j)) + & 
	& 0.25*(x(i+1)-x(i-1))*(uy1(i,j+1)**2-uy1(i,j-1)**2))*dt + uy1(i,j)*Vm(i,j)
      end do
    end do

    call RESOL_THOMAS(Nx,Ny,Dg2,Lw2,Up2,SMx,usx)
    call RESOL_THOMAS(Nx,Ny,Dg2,Lw2,Up2,SMy,usy)

    ! ========================= Etape 2 : projection ============================
    ! ===========================================================================

    ! A l'intérieur
    do i=1,Nx-1
      do j=1,Ny-1
	SMp(i,j) = 0.25*(y(j+1)-y(j-1))*((usx(i,j)+usx(i+1,j))-(usx(i,j)+usx(i-1,j))) + &
	& 0.25*(x(i+1)-x(i-1))*((usy(i,j)+usy(i,j+1))-(usy(i,j)+usy(i,j-1)))
	SMp(i,j) = SMp(i,j)/dt
      end do
    end do	

    ! pour les CLs Neumann hors coins Neumann

    do j=1,Ny-1
      SMp(0,j) = 0.25*(y(j+1)-y(j-1))*usx(1,j)/dt
      SMp(Nx,j) = -0.25*(y(j+1)-y(j-1))*usx(Nx-1,j)/dt
    end do

    do i=1,Nx-1
      SMp(i,0) = 0.25*(x(i+1)-x(i-1))*usy(i,1)/dt
      SMp(i,Ny) = -0.25*(x(i+1)-x(i-1))*usy(i,Ny-1)/dt
    end do

    call RESOL_THOMAS(Nx,Ny,Dg3,Lw3,Up3,SMp,p)

    ! ========================= Etape 3 : résolution directe de u(n+1) ==========
    ! ===========================================================================

    do i=1,Nx-1
      do j=1,Ny-1
	ux2(i,j) = usx(i,j) - 0.25*(dt*(y(j+1)-y(j-1)))*(p(i+1,j)-p(i-1,j))/Vm(i,j)
	uy2(i,j) = usy(i,j) - 0.25*(dt*(x(i+1)-x(i-1)))*(p(i,j+1)-p(i,j-1))/Vm(i,j)
      end do
    end do	

    do i=0,Nx
      ux2(i,0) = 0.
      uy2(i,0) = 0.
      ux2(i,Ny) = U0
      uy2(i,Ny) = 0.	
    end do

    do j=1,Ny-1
      ux2(0,j) = 0.
      uy2(0,j) = 0.
      ux2(Nx,j) = 0.
      uy2(Nx,j) = 0.	
    end do
    
    !Calcul de l'erreur max relative sur l'amplitude la vitesse entre chaque pas de temps

    do i=0,Nx
      do j=0,Ny
	sup(n,1) = max(sup(n,1),abs(sqrt(ux2(i,j)**2+uy2(i,j)**2)-sqrt(ux1(i,j)**2+uy1(i,j)**2)))		
      enddo
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Sauvegarde dans un fichier graçe au module mod_save_vtk.f90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (n > 0) then

      inc = inc+1

      if (inc == Nsave) then

	inc = 0
	inc2 = inc2+1	    	  

	!call stream_compute(ux2,uy2) (ne calcule correctement)

	do i=0,Nx
	  do j=0,Ny

	    !if (Psi(i,j) > sup(n,3)) then
	      !imax = x(i)
	      !jmax = y(j)
	    ! endif

	    sup(n,2) = max(sup(n,2),abs(p(i,j)))
	    !sup(n,3) = max(sup(n,3),abs(Psi(i,j)))	    	

	    enddo
	enddo

	call save_VTK(ech,ux2,uy2,p,x,y,n,Nx,Ny,Psi)	  

	write(*,'(a,x,i7.7,x,f13.6)') 'iteration :',n,real(n*dt)
	write(*,'(a,x,e13.6)') 'Err max vitesse :',sup(n,1)
	write(*,'(a,x,f13.6)') 'max pression :',sup(n,2)
	!write(*,'(a,x,e16.10)') 'max stream function :',sup(n,3)

      end if

    end if

    ux1 = ux2
    uy1 = uy2

  end do
  
  call save_VTK(ech,ux2,uy2,p,x,y,n,Nx,Ny,Psi)
  
  call cpu_time(tps2)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Récapitulatif à la fin de la simulation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,*) '==================================================================================================='
  write(*,'(a,i3.3,x,a,f16.6,a)') "fin simulation Chorin n°",ech,'temps "simulation"', (n-1)*dt,'s'
  !write(*,'(a,f10.6,x,a,f10.6)') 'vortex primaire en x = ',imax,'y = ',jmax
  write(*,'(a,x,i7.7)') 'nombre iteration :',n
  write(*,'(a,x,f16.6)') 'temps de calcul :', tps2-tps1
  write(*,*) 'U0 =', U0
  write(*,*) 'Re = ', Re
  write(*,*) 'dt =', dt,'(s)'
  write(*,*) '==================================================================================================='
  write(*,*)


  open(2,file='Nech.data',form='formatted',position='rewind')
  write(2,'(i3.3)') ech+1
  close(2)

  deallocate(Up2,Dg2,Lw2)
  deallocate(Up3,Dg3,Lw3)
  deallocate(UpPsi,DgPsi,LwPsi)
  deallocate(usx,usy,ux1,ux2,uy1,uy2,p,Psi,SMx,SMy,SMp,SMs)
  deallocate(sup)

end subroutine
