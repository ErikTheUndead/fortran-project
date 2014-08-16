MODULE mod_parameter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!module comprenant toutes les variables communes entre les différentes subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
	
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Variables et paramètres communs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER(kind=4) :: Nx,Ny ! Nombre de point dans les directions x et y
  INTEGER :: maillage_type ! type de maillage 0 = régulier &  1 = irrégulier
  INTEGER :: Eq ! 0 = Poisson, 1 = Helmholtz, 2 = Chorin colocatif et 3 = goda colocatif

  REAL(kind=8) :: Lx,Ly,Dmin ! domaine inf et sup (Lx = Ly dans notre cas), et maille la plus grossière
  REAL(kind=8) :: dt,cv ! pas de  temps et critère de convergence

  REAL(kind=8), allocatable, dimension(:) :: x,y ! définition du maillage	
  REAL(kind=8), allocatable, dimension(:,:) :: Vm ! Volume de la maille (i,j) 
  REAL(kind=8), allocatable, dimension(:,:,:) :: Up,Dg,Lw ! Matrice Upper, Matrice Diag, Matrice Lower

  !!!!!!!!!!!!!!!!!!!!!
  !Poisson et Helmholtz
  !!!!!!!!!!!!!!!!!!!!!

  INTEGER :: CL_type, Neu_D ! Helmholtz.f90 et Poisson.f90 type de CL aux bords du domaine, 0 = Dirichlet et 1 = Neumann, Neu_D = 1 pts Dirichlet / Neumann 4 pts
  REAL(kind=8) :: cx,cy ! cx = Pi/Nx et cy = Pi/Ny
  REAL(KIND=8), parameter :: Pi = acos(-1.)

  REAL(kind=8), allocatable, dimension(:,:) :: f ! Fonction que doit satisfaire les CL de Neumann aux bords  
  REAL(kind=8), allocatable, dimension(:,:) :: Tv,T2 ! Température sol analytique & Température recherché
  REAL(kind=8), allocatable, dimension(:,:) :: T1 ! température au temps n-1 pour Helmholtz
  REAL(kind=8), allocatable, dimension(:,:) :: Sv,S2 ! Source sol analytique & source formulé en VF
  REAL(kind=8), allocatable, dimension(:) :: erreur !

  !!!!!!!!!!!!!!
  !Navier-Stokes
  !!!!!!!!!!!!!!

  INTEGER(kind=4) :: itsave ! indice de l'enregistrement pour save_vtk
  INTEGER(KIND=8) :: iter, Nsave ! nombre d'itération max et intervalle de sauvegarde
  INTEGER :: ech ! numéro de l'echantillon de la donnée pour la sauvegarde des fichiers avec mod_save_vtk.f90

  REAL(kind=8) :: Temps, cfl ! Temps = nombre d'itération x dt & critère CFL
  REAL(KIND=8) :: Re ! Nombre de Reynolds
  REAL(KIND=8) :: U0 ! Vitesse d'entrainement

  REAL(kind=8), allocatable, dimension(:,:) :: usx ! composante sur x de u*
  REAL(kind=8), allocatable, dimension(:,:) :: usy ! composante sur y de u*
  REAL(kind=8), allocatable, dimension(:,:) :: ux1 ! composante sur x au temps(n) de u
  REAL(kind=8), allocatable, dimension(:,:) :: uy1 ! composante sur y au temps(n) de u
  REAL(kind=8), allocatable, dimension(:,:) :: ux2 ! composante sur x au temps(n+1) de u
  REAL(kind=8), allocatable, dimension(:,:) :: uy2 ! composante sur y au temps(n+1) de u
  REAL(kind=8), allocatable, dimension(:,:) :: Psi ! fonction de courant
  REAL(kind=8), allocatable, dimension(:,:) :: p,phi,p1,p2 ! pression au temps (n) (Chorin), pression corrigée (goda), p(n) et p(n+1) (goda)
  REAL(kind=8), allocatable, dimension(:,:) :: SMx,SMy,SMp,SMs ! second membre pour u*|x et v*|y, calcul de la pression et Stream function
  REAL(kind=8), allocatable, dimension(:,:) :: sup ! Calcul de l'erreur à chaque pas de temps

  REAL(kind=8), allocatable, dimension(:,:,:) :: Up2,Dg2,Lw2 ! Matrice Upper, Matrice Diag, Matrice Lower (vitesse pour NS)
  REAL(kind=8), allocatable, dimension(:,:,:) :: Up3,Dg3,Lw3 ! Matrice Upper, Matrice Diag, Matrice Lower (pression pour NS)
  REAL(kind=8), allocatable, dimension(:,:,:) :: UpPsi,DgPsi,LwPsi ! Matrice Upper, Matrice Diag, Matrice Lower (fonction de courant Psi pour NS)

END MODULE
