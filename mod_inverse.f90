!#############################################################################
!#############################################################################
! 
!	TITRE :			     $$$ mod_inverse.f $$$ 
!     -----
!-----------------------------------------------------------------------------
!	AUTEURS : UKAEA, HARWELL						 
!	-------
!-----------------------------------------------------------------------------                    
!	OBJET :  INVERSION DE MATRICE
!     -----	
!		   
!-----------------------------------------------------------------------------
!	CONTENU :
!	-------
!	SUBROUTINE INVERSE
!	
!#############################################################################
!#############################################################################
	
	
	MODULE MOD_INVERSE
	
	contains
	
	
	SUBROUTINE inverse(A)
!......................................................................
!   CETTE SUBROUTINE CALCULE LA MATRICE A-1 INVERSE DE LA MATRICE A.
!   LA MATRICE A EST DETRUITE ET REMPLACEE PAR A-1.
!   ----> N : DIMENSION DE A.
!   ----> A : MATRICE A INVERSER.
!   ----> WORK : TABLEAU DE TRAVAIL.
!......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision, dimension(:,:),intent(INOUT) :: A
      integer  :: M
!      INTEGER  :: IW(1000)
!      double precision, ALLOCATABLE :: 
!     &  WORK(:) !WORK(:,:)
!      ALLOCATE(
!     &         WORK(size(A,1)) !WORK(size(A,1),size(A,1))
!     &        )
	
      M=size(A,1)
      CALL MB01CD(A,M) !,N) !,IW) !,WORK)

!      CALL MB02CD(A)
!      DEALLOCATE(
!     &  WORK
!     &          )

      END SUBROUTINE inverse
!--------------------------------------------------------------------------
!/     ADD NAME=MB01CD          HSL     VAX
!######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
!######ALIAS MB01CD MB01BD MB01AD
!###### CALLS   FM02
      SUBROUTINE MB01CD(A,M) !,IA) !,IND) !,C)
!  STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
!x      DOUBLE PRECISION A,AMAX,C,DIV,STO,W,W1,FM02AD,ZERO,ONE
!      double precision, dimension(:,:),intent(INOUT) :: A
      INTEGER IND,IA
      DOUBLE PRECISION A,AMAX,C,DIV,STO,W,W1,ZERO,ONE
      DIMENSION A(M,M),IND(M),C(M)


!      DIMENSION A(IA,1),IND(1),C(1)
!x      EXTERNAL MB01ED
!x      COMMON/MB01DD/LP,IFLAG
      DATA LP/6/
      DATA ZERO,ONE/0.0D0,1.0D0/
      IA=M
      IFLAG=0
      IF(M-1)50,2,3
2     IF(A(1,1).EQ.ZERO)GO TO 60
      A(1,1)=ONE/A(1,1)
      GO TO 99
    3 M1=M-1
      AMAX=ZERO
      DO 32 I=1,M
      IND(I)=I
      IF(DABS(A(I,1))-DABS(AMAX))32,32,31
   31 AMAX=A(I,1)
      IMAX=I
   32 CONTINUE
      IF(AMAX.EQ.ZERO)GO TO 60
      DO 41 J=1,M1
      IF(IMAX-J)35,35,33
   33 IW=IND(IMAX)
      IND(IMAX)=IND(J)
      IND(J)=IW
      DO 34 K=1,M
      W=A(IMAX,K)
      A(IMAX,K)=A(J,K)
      A(J,K)=W
   34 CONTINUE
   35 J1=J+1
      IF(J.EQ.1)GO TO 38
   36 DO 37 I=J1,M
      A(J,I)=A(J,I)-FM02AD(J-1,A(J,1),IA,A(1,I),1)
   37 CONTINUE
   38 DIV=AMAX
      AMAX=ZERO
      DO 40 I=J1,M
      A(I,J)=A(I,J)/DIV
      A(I,J+1)=A(I,J+1)-FM02AD(J,A(I,1),IA,A(1,J+1),1)
      IF(DABS(A(I,J1))-DABS(AMAX))40,40,39
   39 AMAX=A(I,J1)
      IMAX=I
   40 CONTINUE
      IF(AMAX.EQ.ZERO)GO TO 60
   41 CONTINUE
      DO 13 I1=1,M1
      I=M+1-I1
      I2=I-1
      DO 11 J1=1,I2
      J=I2+1-J1
      J2=J+1
      W1=-A(I,J)
      IF(I2-J2)10,9,9
    9 W1=W1-FM02AD(I2-J2+1,A(J2,J),1,C(J2),1)
   10 C(J)=W1
   11 CONTINUE
      DO 12 K=1,I2
      A(I,K)=C(K)
   12 CONTINUE
   13 CONTINUE
      DO 22 I1=1,M
      I=M+1-I1
      I2=I+1
      W=A(I,I)
      DO 20 J=1,M
      IF(I-J)14,15,16
   14 W1=ZERO
      GO TO 17
   15 W1=ONE
      GO TO 17
   16 W1=A(I,J)
   17 IF(I1-1)19,19,18
   18 W1=W1-FM02AD(M-I2+1,A(I,I2),IA,A(I2,J),1)
   19 C(J)=W1
   20 CONTINUE
      DO 21 J=1,M
      A(I,J)=C(J)/W
   21 CONTINUE
   22 CONTINUE
      DO 26 I=1,M
   23 IF(IND(I)-I)24,26,24
   24 J=IND(I)
      DO 25 K=1,M
      STO=A(K,I)
      A(K,I)=A(K,J)
      A(K,J)=STO
   25 CONTINUE
      ISTO=IND(J)
      IND(J)=J
      IND(I)=ISTO
      GO TO 23
   26 CONTINUE
      GO TO 99
50    IF(LP.GT.0)WRITE(LP,55)
55    FORMAT(51H ERROR RETURN FROM MB01CD BECAUSE M IS NOT POSITIVE)
      IFLAG=1
      GO TO 99
60    IF(LP.GT.0)WRITE(LP,65)
65    FORMAT(52H ERROR RETURN FROM MB01CD BECAUSE MATRIX IS SINGULAR)
      IFLAG=2
   99 RETURN
      END SUBROUTINE MB01CD
!--------------------------------------------------------------------------      
      DOUBLE PRECISION FUNCTION FM02AD(N,A,IA,B,IB)
      DOUBLE PRECISION R1,A,B
      DIMENSION A(*),B(*)

!    N   THE LENGTH OF THE VECTORS (IF N<= 0  FM02AD = 0)
!    A   THE FIRST VECTOR
!    IA  SUBSCRIPT DISPLACEMENT BETWEEN ELEMENTS OF A
!    B   THE SECOND VECTOR
!    IB  SUBSCRIPT DISPLACEMENT BETWEEN ELEMENTS OF B
!    FM02AD  THE RESULT

      R1=0D0
      IF(N.LE.0) GO TO 2
      JA=1
      IF(IA.LT.0) JA=1-(N-1)*IA
      JB=1
      IF(IB.LT.0) JB=1-(N-1)*IB
      I=0
    1 I=I+1
      R1=R1+A(JA)*B(JB)
      JA=JA+IA
      JB=JB+IB
      IF(I.LT.N) GO TO 1
    2 FM02AD=R1
!      RETURN
      END FUNCTION FM02AD
      

!_______________________________________________________________________      


	END MODULE MOD_INVERSE


     
