!############################################################################
!############################################################################
!  
!	 TITRE :               $$$ mod_thomas.f $$$
!      -----
!-----------------------------------------------------------------------------
!	AUTEUR : V.DELGADO - A.REDONDO						 
!	------
!-----------------------------------------------------------------------------
!	OBJET :	 MODULE DE RESOLUTION DIRECTE - ALGORITHME DE THOMAS
!	-----				
!
!#############################################################################
!#############################################################################

	MODULE MOD_THOMAS	
	USE MOD_INVERSE
	CONTAINS


!************************************************************************
!************************************************************************
!======================SUBROUTINE INIT_THOMAS=============================
!************************************************************************
!************************************************************************
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  ENTREE : ------------------------------
!          
!	
!
!------------------------------- SORTIE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>		
!     
!	
!   
!------------------------------------------------------------------------

	SUBROUTINE INIT_THOMAS(N,M,AD,AL,AU)
	
	
	!------------ DECLARATION DES VARIABLES -----------------------

	implicit double precision (A-H,O-Z)
	double precision :: AD
	
	

	!VARIABLES GLOBALES:

      dimension AD(0:N,0:N,0:M),AL(0:N,0:N,0:M),AU(0:N,0:N,0:M),&
     &Mk(0:N,0:N,0:M)    	

	
	!--------------------------------------------------------------
	
	!Inversion des opérateurs situés sur la diagonale
	
	AL(:,:,0)=AD(:,:,0)
	call inverse(AL(:,:,0))	
	AU(:,:,0)=matmul(AL(:,:,0),AU(:,:,0))
	
	do k=1,M
		!Modification de Dk->D'k		
		AD(:,:,k)=AD(:,:,k)-matmul(AL(:,:,k),AU(:,:,k-1))
		call inverse(AD(:,:,k))

	    !Calcul de Mk		
		AU(:,:,k)=matmul(AD(:,:,k),AU(:,:,k))
	enddo  	
	

	END SUBROUTINE INIT_THOMAS

!************************************************************************
!************************************************************************
!======================SUBROUTINE RESOL_THOMAS=============================
!************************************************************************
!************************************************************************
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  ENTREE : ------------------------------
!          
!	
!
!------------------------------- SORTIE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>		
!     
!	
!   
!------------------------------------------------------------------------

	SUBROUTINE RESOL_THOMAS(N,M,AD,AL,AU,MAT_S,SOL)
	
	
	!------------ DECLARATION DES VARIABLES -----------------------

	implicit double precision (A-H,O-Z)	

	!VARIABLES GLOBALES:
	
	double precision, dimension (:,:), allocatable ::MAT_Sb
	double precision,dimension(:,:),allocatable::Beta_k,WORK

	double precision :: MAT_S,SOL

      dimension AD(0:N,0:N,0:M),AL(0:N,0:N,0:M)
	dimension AU(0:N,0:N,0:M),Mk(0:N,0:N,0:M)
	dimension MAT_S(0:N,0:M),SOL(0:N,0:M)
	allocate(MAT_Sb(0:N,0:M))
	allocate(Beta_k(0:N,0:M))
	allocate(WORK(0:N,0:M))

	! Ecriture du second membre modifié
	
	MAT_Sb=MAT_S

	!j=0
	Beta_k(:,0)=matmul(AL(:,:,0),MAT_S(:,0))
		

	do j=1,M
	WORK(:,j)=MAT_Sb(:,j)-matmul(AL(:,:,j),Beta_k(:,j-1))
	Beta_k(:,j)=matmul(AD(:,:,j),WORK(:,j))		
	enddo

	!Résolution :

	SOL=0.d0

	!j=M	
	SOL(:,M)=Beta_k(:,M)
	
	
	do j=M-1,0,-1		
	SOL(:,j)=Beta_k(:,j)-matmul(AU(:,:,j),SOL(:,j+1))	
	enddo


	DEALLOCATE(MAT_Sb,Beta_k,WORK)
	
	END SUBROUTINE RESOL_THOMAS   
 

	END MODULE MOD_THOMAS
