SUBROUTINE init(temp, psi, vort, u, v, Nstep_old)
USE Parametres

implicit none
    
double precision ,dimension(Nx1+Nx2+1,Ny1+Ny2+1)  :: temp , psi ,vort , u, v     
double precision :: tmp    
integer :: Nstep_old
logical :: log1
    
integer:: ix,iy
!
!
!
!    _______________________________________
!    |                                      |
!    |                                      |
!    |       t=0         U=1                |
!    |                                      |
!    |__________    _  _  _  _  _  _  _  _  |
!               |_________U=0_______________|
!                          
!
!   CONDITION INITIALES
!=======================================
if(isuite.eq.0) then  
!---------------------- ! Pas de suite de calcul
  Nstep_old=0          ! On initialise les variables normalement 
  temp=0.d0
  vort=0.d0
  v=0.d0
  
  ! Condition initiale sur la premiere composante de vitesse
  do iy=1,Ny1+Ny2+1
    if(iy.le.Ny1) then
      do ix=1,Nx2+1
        u(ix+Nx1,iy)=0.d0
      enddo
    else
      do ix=1,Nx1+Nx2+1
        u(ix,iy)=1.d0
      enddo
    endif 
  enddo
  
  ! On deduit psi a partir de u
  
  do iy=1,Ny1+Ny2+1
    if(iy.le.Ny1) then
      do ix=1,Nx2+1
        psi(ix+Nx1,iy)=0.d0
      enddo
    else
      do ix=1,Nx1+Nx2+1
        psi(ix,iy)=(iy-Ny1-1)*dy
      enddo
    endif 
  enddo
!=======
else                !On recupere le fichier res_resulat de l'ancien calcul pour lancer 
!=======             ! une suite de calcul

  inquire (file="res_suite",exist=log1)
  ! Lecture du fichier suite
  if (log1.eqv..true.) then
    open(unit=Nstep+501, file='res_suite',status='old',&
               access='sequential', form='formatted', action='read')
    read(Nstep+501,*) Nstep_old
    if(Nstep_old.ge.Nstep)then
       write(*,*) 'L option  suite de calcul est activee dans parametres  (isuite =1)'
       write(*,*) 'Le nouveau nb d iteration doit etre superieur a l ancien nb d iteration'
       write(*,*) 'l ancien nb d iteration est :  ', Nstep_old, ' et vous avez entré :' , Nstep , '!!!'
       stop
    endif
    do ix=1,Nx1+Nx2+1
       do iy=1,Ny1+Ny2+1
         read(Nstep+501,*) tmp , tmp, u(ix, iy), v(ix, iy), psi(ix, iy), vort(ix, iy), temp(ix, iy) 
       enddo
       read(Nstep+501,*) 
    enddo
    close(501+Nstep)
  else
       write(*,*) 'L option  suite de calcul est activee dans parametres (isuite =1)'
    write(*,*) 'Attention !! Pas de fichier pour lancer une suite de calcul' 
    stop
  endif
endif  
    
END SUBROUTINE init
