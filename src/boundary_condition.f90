SUBROUTINE bc(temp, psi, vort, u, v)
USE Parametres

implicit none
    
double precision ,dimension(Nx1+Nx2+1,Ny1+Ny2+1)  :: temp , psi ,vort , u, v     
     
integer:: ix,iy
!
!
!
!    _______________________________________
!    |                  B5                  |
!    |                                      |
!  B1|                                      |B6
!    |                                      |
!    |__________ B3                         |
!       B2      |___________________________|
!                          B4
!
!   CONDITION LIMITE SUR LA PAROI B1
!=======================================
    do iy=1,Ny2+1
       temp(1,iy+Ny1)=0.D0
       psi(1,iy+Ny1)=dy*(iy-1)
       vort(1,iy+Ny1)=0.d0
       u(1,iy+Ny1)=1.d0
       v(1,iy+Ny1)=0.d0
    enddo

!   CONDITION LIMITE SUR LA PAROI B2
!=======================================
    do ix=1,Nx1+1
       temp(ix,Ny1+1)=0.d0
       psi(ix,Ny1+1 )=0.d0
       vort(ix, Ny1+1 )=2.d0*(psi(ix,(Ny1+1)+1)-dy)/dy**2
       u(ix,Ny1+1 )=0.d0
       v(ix,Ny1+1 )=0.d0
    enddo


!   CONDITION LIMITE SUR LA PAROI B3
!=======================================
    do iy=1,Ny1
       temp(Nx1+1,iy)=1.d0
       psi(Nx1+1,iy )=0.d0
       vort(Nx1+1,iy)=2.d0*psi((Nx1+1)+1,iy)/dx**2
       u(Nx1+1,iy )=0.d0
       v(Nx1+1,iy )=0.d0
    enddo


!   CONDITION LIMITE SUR LA PAROI B4
!=======================================
    do ix=1,Nx2
       temp(ix+Nx1+1,1)=0.d0
       psi(ix+Nx1+1,1)=0.d0
       vort(ix+Nx1+1,1)=2*psi(ix+Nx1+1,2)/dy**2
       u(ix+Nx1+1,1)=0.d0
       v(ix+Nx1+1,1)=0.d0
    enddo


!   CONDITION LIMITE SUR LA PAROI B5
!======================================
    do ix=1, Nx1+Nx2
       temp(ix+1,Ny1+Ny2+1)=0.d0
       psi(ix+1,Ny1+Ny2+1)=dy*Ny2
       vort(ix+1,Ny1+Ny2+1)=2.d0*(psi(ix+1,(Ny1+Ny2+1)-1)-psi(ix+1,Ny1+Ny2+1)+dy)/dy**2
       u(ix+1,Ny1+Ny2+1)=1.d0
       v(ix+1,Ny1+Ny2+1)=0.d0
    enddo


!   CONDITION LIMITE SUR LA PAROI B6
!=======================================
    do iy=1,Ny1+Ny2-1
       temp(Nx1+Nx2+1,iy+1)=temp((Nx1+Nx2+1)-1,iy+1)
       psi(Nx1+Nx2+1,iy+1 )=psi((Nx1+Nx2+1)-1,iy+1)
       vort(Nx1+Nx2+1,iy+1 )=vort((Nx1+Nx2+1)-1,iy+1)
       u(Nx1+Nx2+1,iy+1 )=u((Nx1+Nx2+1)-1,iy+1)
       v(Nx1+Nx2+1,iy+1 )=v((Nx1+Nx2+1)-1,iy+1)
    enddo


END SUBROUTINE bc
