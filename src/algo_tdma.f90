!==============================================
! Resolution d'un systeme lineaire tridiagonal
! algorithme TDMA ou thoma
! N la taille du systeme
! F solution cherchee
! A, B, C les diag
! MF=D    
!    (a1 b1 ...0...0  )              
!    (c2 a2 .      .  ) 
! M =(0  .     .   0  )
!    (.     .     bN-1)    
!    (0....0...cN aN  )
!==============================================
!***********************************************
!         TDMA Algorithm 
!***********************************************
 
subroutine TDMA(A0,B0,C0,D0,F0,N)

implicit none 
integer :: k, N
double precision  :: A0(N),B0(N),C0(N),D0(N),F0(N)
double precision , dimension (:), allocatable:: P, Q 

allocate (P(N))
allocate (Q(N))
F0=0.d0
!Calcul de P et Q
C0(1)=0.d0
B0(N)=0.d0
P(1)=B0(1)/A0(1)
Q(1)=D0(1)/A0(1)

do k = 2,N
 P(k)=B0(k)/(A0(k)-C0(k)*P(k-1))
 Q(k)=(D0(k)+C0(k)*Q(k-1))/(A0(k)-C0(k)*P(k-1))
enddo
F0(N)=Q(N)
do k=N-1,1,-1
 F0(k)=P(k)*F0(k+1)+Q(k)
enddo

deallocate(P,Q)

end
!===============================================
