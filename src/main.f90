!=====================================
! Re:  Reynolds number
! Ri:  Richardson number
! Pr:  Prandtl number
! Er:  Expansion rate (Er-1: backstep height)
! Xs:  Position on the step from the inlet
! Xt:  Distance between the inlet and the outlet
! dt:  time step
! dx:  X step 
! dy:  Y step 
! Nt:  Number of time iterations
! Nx:  Number of cells in the X direction
! Ny:  Number of cells in the Y direction
! Fox:  Fourier number following x (a*dt/dx)
! Foy:  Fourier number following y (a*dt/dy)
! Cx:  Courant number following x (a*dt/dx)
! Cy:  Courant number following y (a*dt/dy)
! temp :  Temperature T(i,j)
! vort:  Vorticity Vc(i,j)
! psi: Courant function Phi(i,j)
! 
!====================================


program backstep 

use Parametres

implicit none 

integer :: Nt, ixs, iys
!double precision , dimension (:,:) :: T, Vc, Phi
!coordonnee des centres des cellules
double precision , dimension (Nx1+Nx2+1):: xx
double precision , dimension (Ny1+Ny2+1):: yy
!champs de vitesses
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: u, v
!Temperature
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: temp
!fonction de courant
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: psi
!copie de la fonction de courant
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: copyPsi
!vorticite
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: vort
!tableaux a utiliser pour tdma
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: atemp, btemp, ctemp, dtemp  ! temperature
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: avort, bvort, cvort, dvort  ! vorticite 
double precision, dimension(Nx1+Nx2+1,Ny1+Ny2+1):: apsi, bpsi, cpsi, dpsi  ! psi 

double precision :: diff
integer :: ix ,iy,  istep, comp, Nstep_old
character (len=4) :: base
character (len=80) :: int2char, filename

!*** Remplissage du tableau de la position
!*****************************************
do ix=1,Nx1+Nx2+1
  xx(ix)=(ix-1)*dx
enddo
do iy=1,Ny1+Ny2+1
  yy(iy)=(iy-1)*dy
enddo

!***Coordonnée des sondes 
!==========================
do ix=1,Nx1+Nx2+1
  if(abs(xx(ix)-Xsonde).le.0.5d0*dx) then
    do iy=1,Ny1+Ny2+1
      if(abs(yy(iy)-Ysonde).le.0.5d0*dy) then
         ixs=ix
         iys=iy
      endif
    enddo
  endif
enddo
!=============
if(.true.) then
!========INITIALISATION===============
  call init(temp,psi,vort, u,v, Nstep_old)
!==================REMPLISSAGE DES TABLEAUX DE LA VITESSE============
!--------------- 
open (unit=11, file='sonde.dat')
write(11,*) '# 1- u  2- v  3- psi  4-omega  5- temperature'
do istep=Nstep_old+1, Nstep
  call bc(temp, psi, vort, u, v)  ! Les conditions aux limites
  do iy=2,Ny1+Ny2
   if (iy.le.Ny1+1) then
     do ix=2,Nx2     
      u(ix+Nx1,iy)=(psi(ix+Nx1,iy+1)-psi(ix+Nx1,iy-1))/(2.d0*dy)
      v(ix+Nx1,iy)=-(psi(ix+Nx1+1,iy)-psi(ix+Nx1-1,iy))/(2.d0*dx)                    
     enddo
   else
     do ix=2,Nx1+Nx2     
      u(ix,iy)=(psi(ix,iy+1)-psi(ix,iy-1))/(2.d0*dy)
      v(ix,iy)=-(psi(ix+1,iy)-psi(ix-1,iy))/(2.d0*dx)                    
     enddo
   endif
  enddo
!===============RESOLUTION DE LA TEMPERATURE PAR LA METHODE ADI==========
!=====Remplissage des tableaux a, b c, d de TDMA pour la premiere etape de l'ADI===
 do iy=2,Ny1+Ny2
   if (iy.le.Ny1+1) then
     do ix=2,Nx2     
        atemp(ix+Nx1,iy)=1.d0+Fox       
        btemp(ix+Nx1,iy)=-u(ix+Nx1,iy)*Cx/4.d0 + Fox/2.d0      
        ctemp(ix+Nx1,iy)=u(ix+Nx1,iy)*Cx/4.d0 +Fox/2.d0       
       if(ix.eq.2) then
         dtemp(ix+Nx1,iy)=temp(ix+Nx1,iy)-(Cy*v(ix+Nx1,iy)/4.d0)*&
                          (temp(ix+Nx1,iy+1)-temp(ix+Nx1, iy-1))+&
                          (Foy/2.d0)*(temp(ix+Nx1,iy+1)+temp(ix+Nx1,iy-1)-2.d0*temp(ix+Nx1,iy))+&
                          ctemp(ix+Nx1,iy)*temp(1+Nx1,iy)
 
        elseif(ix.eq.Nx2) then
          dtemp(ix+Nx1,iy)=temp(ix+Nx1,iy)-(Cy*v(ix+Nx1,iy)/4.d0)*&
                          (temp(ix+Nx1,iy+1)-temp(ix+Nx1, iy-1))+&
                          (Foy/2.d0)*(temp(ix+Nx1,iy+1)+temp(ix+Nx1,iy-1)-2.d0*temp(ix+Nx1,iy))+&
                          btemp(ix+Nx1,iy)*temp(1+Nx1+Nx2,iy)

        else
          dtemp(ix+Nx1,iy)=temp(ix+Nx1,iy)-(Cy*v(ix+Nx1,iy)/4.d0)*&
                          (temp(ix+Nx1,iy+1)-temp(ix+Nx1, iy-1))+&
                          (Foy/2.d0)*(temp(ix+Nx1,iy+1)+temp(ix+Nx1,iy-1)-2.d0*temp(ix+Nx1,iy))
        
        endif
     enddo
   else
     do ix=2,Nx1+Nx2     
        atemp(ix,iy)=1.d0+Fox       
        btemp(ix,iy)=-u(ix,iy)*Cx/4.d0 + Fox/2.d0      
        ctemp(ix,iy)=u(ix,iy)*Cx/4.d0 +Fox/2.d0       
        if(ix.eq.2) then
           dtemp(ix,iy)=temp(ix,iy)-(Cy*v(ix,iy)/4.d0)*&
                          (temp(ix,iy+1)-temp(ix, iy-1))+&
                          (Foy/2.d0)*(temp(ix,iy+1)+temp(ix, iy-1)-2.d0*temp(ix,iy))+&
                          ctemp(ix,iy)*temp(1,iy)
         elseif(ix.eq.Nx1+Nx2) then
           dtemp(ix,iy)=temp(ix,iy)-(Cy*v(ix,iy)/4.d0)*&
                          (temp(ix,iy+1)-temp(ix, iy-1))+&
                          (Foy/2.d0)*(temp(ix,iy+1)+temp(ix, iy-1)-2.d0*temp(ix,iy))+&
                          btemp(ix,iy)*temp(Nx1+Nx2+1,iy)
 
         else
           dtemp(ix,iy)=temp(ix,iy)-(Cy*v(ix,iy)/4.d0)*&
                          (temp(ix,iy+1)-temp(ix, iy-1))+&
                          (Foy/2.d0)*(temp(ix,iy+1)+temp(ix, iy-1)-2.d0*temp(ix,iy))
 
         endif


     enddo
   endif
 enddo
!=========Resolution de T(n+1/2)==
 do iy=2,Ny1+Ny2
   if (iy.le.Ny1+1) then
     call TDMA(atemp((Nx1+2):(Nx1+Nx1),iy), btemp((Nx1+2):(Nx1+Nx2),iy),&
               ctemp((Nx1+2):(Nx1+Nx2),iy), dtemp((Nx1+2):(Nx1+Nx2),iy),&
               temp((Nx1+2):(Nx1+Nx2),iy),Nx2-1)                      
   else
      call TDMA(atemp(2:(Nx1+Nx2),iy), btemp(2:(Nx1+Nx2),iy),&
               ctemp(2:(Nx1+Nx2),iy), dtemp(2:(Nx1+Nx2),iy),&
               temp(2:(Nx1+Nx2),iy),Nx1+Nx2-1)                      
   endif
 enddo

!====Remplissage des tableaux a, b c, d de TDMA pour la deuxieme  etape de l'ADI==
 do ix=2,Nx1+Nx2
   if (ix.le.Nx1+1) then
     do iy=2,Ny2     
        atemp(ix,iy+Ny1)=1.d0+Foy       
        btemp(ix,iy+Ny1)=-v(ix,iy+Ny1)*Cy/4.d0 + Foy/2.d0      
        ctemp(ix,iy+Ny1)=v(ix,iy+Ny1)*Cy/4.d0 +Foy/2.d0       
        if(iy.eq.2) then
           dtemp(ix,iy+Ny1)=temp(ix,iy+Ny1)-(Cx*u(ix,iy+Ny1)/4.d0)*&
                          (temp(ix+1,iy+Ny1)-temp(ix-1, iy+Ny1))+&
                          (Fox/2.d0)*(temp(ix+1,iy+Ny1)+temp(ix-1, iy+Ny1)-2.d0*temp(ix,iy+Ny1))+&
                           ctemp(ix,iy+Ny1)*temp(ix, Ny1+1)
        elseif(iy.eq.Ny2) then
            dtemp(ix,iy+Ny1)=temp(ix,iy+Ny1)-(Cx*u(ix,iy+Ny1)/4.d0)*&
                          (temp(ix+1,iy+Ny1)-temp(ix-1, iy+Ny1))+&
                          (Fox/2.d0)*(temp(ix+1,iy+Ny1)+temp(ix-1, iy+Ny1)-2.d0*temp(ix,iy+Ny1))+&
                           btemp(ix,iy+Ny1)*temp(ix, Ny1+Ny2+1)
 
        else
           dtemp(ix,iy+Ny1)=temp(ix,iy+Ny1)-(Cx*u(ix,iy+Ny1)/4.d0)*&
                          (temp(ix+1,iy+Ny1)-temp(ix-1, iy+Ny1))+&
                          (Fox/2.d0)*(temp(ix+1,iy+Ny1)+temp(ix-1, iy+Ny1)-2.d0*temp(ix,iy+Ny1))
        endif
     enddo
   else
     do iy=2,Ny1+Ny2     
        atemp(ix,iy)=1.d0+Foy       
        btemp(ix,iy)=-v(ix,iy)*Cy/4.d0 + Foy/2.d0      
        ctemp(ix,iy)=v(ix,iy)*Cy/4.d0 +Foy/2.d0       
        if(iy.eq.2) then
            dtemp(ix,iy)=temp(ix,iy)-(Cx*u(ix,iy)/4.d0)*&
                          (temp(ix+1,iy)-temp(ix-1, iy))+&
                          (Fox/2.d0)*(temp(ix+1,iy)+temp(ix-1, iy)-2.d0*temp(ix,iy))+&
                           ctemp(ix,iy)*temp(ix,1)
         elseif(iy.eq.Ny1+Ny2) then
            dtemp(ix,iy)=temp(ix,iy)-(Cx*u(ix,iy)/4.d0)*&
                          (temp(ix+1,iy)-temp(ix-1, iy))+&
                          (Fox/2.d0)*(temp(ix+1,iy)+temp(ix-1, iy)-2.d0*temp(ix,iy))+&
                           btemp(ix,iy)*temp(ix,Ny1+Ny2+1)
 
         else
            dtemp(ix,iy)=temp(ix,iy)-(Cx*u(ix,iy)/4.d0)*&
                          (temp(ix+1,iy)-temp(ix-1, iy))+&
                          (Fox/2.d0)*(temp(ix+1,iy)+temp(ix-1, iy)-2.d0*temp(ix,iy))

         endif
     enddo
   endif
 enddo
!========Resolution de la temperature T(n+1)=======
 do ix=2,Nx1+Nx2
   if (ix.le.Nx1+1) then
      call TDMA(atemp(ix,(Ny1+2):(Ny1+Ny2)), btemp(ix,(Ny1+2):(Ny1+Ny2)),&
               ctemp(ix,(Ny1+2):(Ny1+Ny2)), dtemp(ix,(Ny1+2):(Ny1+Ny2)),&
               temp(ix,(Ny1+2):(Ny1+Ny2)),Ny2-1)                      
   else
       call TDMA(atemp(ix,2:(Ny1+Ny2)), btemp(ix,2:(Ny1+Ny2)),&
               ctemp(ix,2:(Ny1+Ny2)), dtemp(ix,2:(Ny1+Ny2)),&
               temp(ix,2:(Ny1+Ny2)),Ny1+Ny2-1)                      
   endif
 enddo

!=====Remplissage des tableaux a, b c, d de TDMA pour la premiere etape de l'ADI===
 do iy=2,Ny1+Ny2
   if (iy.le.Ny1+1) then
     do ix=2,Nx2     
       avort(ix+Nx1,iy)=1.d0+Fox*Pr       
       bvort(ix+Nx1,iy)=-u(ix+Nx1,iy)*Cx/4.d0 + Fox*Pr/2.d0      
       cvort(ix+Nx1,iy)=u(ix+Nx1,iy)*Cx/4.d0 +Fox*Pr/2.d0       
       if(ix.eq.2) then
         dvort(ix+Nx1,iy)=vort(ix+Nx1,iy)-(Cy*v(ix+Nx1,iy)/4.d0)*&
                      (vort(ix+Nx1,iy+1)-vort(ix+Nx1, iy-1))+&
                      (Foy*Pr/2.d0)*(vort(ix+Nx1,iy+1)+vort(ix+Nx1, iy-1)-2.d0*vort(ix+Nx1,iy))&
                      -(Cx*Ri/4.d0)*(temp(ix+Nx1+1,iy)-temp(ix+Nx1-1,iy))+&
                       cvort(ix+Nx1,iy)*vort(1+Nx1,iy)          
 
    
       elseif(ix.eq.Nx2) then
         dvort(ix+Nx1,iy)=vort(ix+Nx1,iy)-(Cy*v(ix+Nx1,iy)/4.d0)*&
                      (vort(ix+Nx1,iy+1)-vort(ix+Nx1, iy-1))+&
                      (Foy*Pr/2.d0)*(temp(ix+Nx1,iy+1)+vort(ix+Nx1, iy-1)-2.d0*vort(ix+Nx1,iy))&
                     -(Cx*Ri/4.d0)*(temp(ix+Nx1+1,iy)-temp(ix+Nx1-1,iy))+&
                      bvort(ix+Nx1,iy)*vort(1+Nx1+Nx2,iy)          
 
      
       else
         dvort(ix+Nx1,iy)=vort(ix+Nx1,iy)-(Cy*v(ix+Nx1,iy)/4.d0)*&
                      (vort(ix+Nx1,iy+1)-vort(ix+Nx1, iy-1))+&
                      (Foy*Pr/2.d0)*(vort(ix+Nx1,iy+1)+vort(ix+Nx1, iy-1)-2.d0*vort(ix+Nx1,iy))&
                     -(Cx*Ri/4.d0)*(temp(ix+Nx1+1,iy)-temp(ix+Nx1-1,iy))
 
        endif
     enddo                  
   else
     do ix=2,Nx1+Nx2     
        avort(ix,iy)=1.d0+Fox*Pr       
        bvort(ix,iy)=-u(ix,iy)*Cx/4.d0 + Fox*Pr/2.d0      
        cvort(ix,iy)=u(ix,iy)*Cx/4.d0 +Fox*Pr/2.d0       
        if(ix.eq.2) then
           dvort(ix,iy)=vort(ix,iy)-(Cy*v(ix,iy)/4.d0)*&
                    (vort(ix,iy+1)-vort(ix, iy-1))+&
                    (Foy*Pr/2.d0)*(vort(ix,iy+1)+vort(ix, iy-1)-2.d0*vort(ix,iy))&
                   -(Cx*Ri/4.d0)*(temp(ix+1,iy)-temp(ix-1,iy))+&
                    ctemp(ix,iy)*vort(1,iy)
        elseif(ix.eq.Nx1+Nx2) then
           dvort(ix,iy)=vort(ix,iy)-(Cy*v(ix,iy)/4.d0)*&
                    (vort(ix,iy+1)-vort(ix, iy-1))+&
                    (Foy*Pr/2.d0)*(vort(ix,iy+1)+vort(ix, iy-1)-2.d0*vort(ix,iy))&
                   -(Cx*Ri/4.d0)*(temp(ix+1,iy)-temp(ix-1,iy))+&
                    bvort(ix,iy)*vort(Nx1+Nx2+1,iy)
        else
           dvort(ix,iy)=vort(ix,iy)-(Cy*v(ix,iy)/4.d0)*&
                    (vort(ix,iy+1)-vort(ix, iy-1))+&
                    (Foy*Pr/2.d0)*(vort(ix,iy+1)+vort(ix, iy-1)-2.d0*vort(ix,iy))&
                    -(Cx*Ri/4.d0)*(temp(ix+1,iy)-temp(ix-1,iy))

        endif

     enddo
   endif
 enddo
!=========Resolution de vort(n+1/2)==
 do iy=2,Ny1+Ny2
   if (iy.le.Ny1+1) then
     call TDMA(avort((Nx1+2):(Nx1+Nx2),iy), bvort((Nx1+2):(Nx1+Nx2),iy),&
               cvort((Nx1+2):(Nx1+Nx2),iy), dvort((Nx1+2):(Nx1+Nx2),iy),&
               vort((Nx1+2):(Nx1+Nx2),iy),Nx2-1)                      
   else
      call TDMA(avort(2:(Nx1+Nx2),iy), bvort(2:(Nx1+Nx2),iy),&
               cvort(2:(Nx1+Nx2),iy), dvort(2:(Nx1+Nx2),iy),&
               vort(2:(Nx1+Nx2),iy),Nx1+Nx2-1)                      
   endif
 enddo

!===Remplissage des tableaux a, b c, d de TDMA pour la deuxieme  etape de l'ADI==
 do ix=2,Nx1+Nx2
   if (ix.le.Nx1+1) then
     do iy=2,Ny2     
        avort(ix,iy+Ny1)=1.d0+Foy*Pr       
        bvort(ix,iy+Ny1)=-v(ix,iy+Ny1)*Cy/4.d0 + Foy*Pr/2.d0      
        cvort(ix,iy+Ny1)=v(ix,iy+Ny1)*Cy/4.d0 +Foy*Pr/2.d0       
        if(iy.eq.2) then
           dvort(ix,iy+Ny1)=vort(ix,iy+Ny1)-(Cx*u(ix,iy+Ny1)/4.d0)*&
                       (vort(ix+1,iy+Ny1)-vort(ix-1, iy+Ny1))+&
                       (Fox*Pr/2.d0)*(vort(ix+1,iy+Ny1)+vort(ix-1, iy+Ny1)-2.d0*vort(ix,iy+Ny1))&
                       -(Cx*Ri/4.d0)*(temp(ix+1,iy+Ny1)-temp(ix-1,iy+Ny1))+&
                        cvort(ix,iy+Ny1)*vort(ix,1+Ny1) 
        elseif(iy.eq.Ny2) then
           dvort(ix,iy+Ny1)=vort(ix,iy+Ny1)-(Cx*u(ix,iy+Ny1)/4.d0)*&
                       (vort(ix+1,iy+Ny1)-vort(ix-1, iy+Ny1))+&
                       (Fox*Pr/2.d0)*(vort(ix+1,iy+Ny1)+vort(ix-1, iy+Ny1)-2.d0*vort(ix,iy+Ny1))&
                       -(Cx*Ri/4.d0)*(temp(ix+1,iy+Ny1)-temp(ix-1,iy+Ny1))+&
                       bvort(ix,iy+Ny1)*vort(ix,Ny1+Ny2+1) 
        else
            dvort(ix,iy+Ny1)=vort(ix,iy+Ny1)-(Cx*u(ix,iy+Ny1)/4.d0)*&
                       (vort(ix+1,iy+Ny1)-vort(ix-1, iy+Ny1))+&
                       (Fox*Pr/2.d0)*(vort(ix+1,iy+Ny1)+vort(ix-1, iy+Ny1)-2.d0*vort(ix,iy+Ny1))&
                       -(Cx*Ri/4.d0)*(temp(ix+1,iy+Ny1)-temp(ix-1,iy+Ny1))
        
        endif
     enddo
   else
     do iy=2,Ny1+Ny2     
        avort(ix,iy)=1.d0+Foy*Pr       
        bvort(ix,iy)=-v(ix,iy)*Cy/4.d0 + Foy*Pr/2.d0      
        cvort(ix,iy)=v(ix,iy)*Cy/4.d0 +Foy*Pr/2.d0       
       if(iy.eq.2) then
        dvort(ix,iy)=vort(ix,iy)-(Cx*u(ix,iy)/4.d0)*&
                     (vort(ix+1,iy)-vort(ix-1, iy))+&
                     (Fox*Pr/2.d0)*(vort(ix+1,iy)+vort(ix-1, iy)-2.d0*vort(ix,iy))&
                    -(Cx*Ri/4.d0)*(temp(ix+1,iy)-temp(ix-1,iy))+&
                     cvort(ix,iy)*vort(ix,1)
       elseif(iy.eq.Ny1+Ny2) then
        dvort(ix,iy)=vort(ix,iy)-(Cx*u(ix,iy)/4.d0)*&
                     (vort(ix+1,iy)-vort(ix-1, iy))+&
                     (Fox*Pr/2.d0)*(vort(ix+1,iy)+vort(ix-1, iy)-2.d0*vort(ix,iy))&
                    -(Cx*Ri/4.d0)*(temp(ix+1,iy)-temp(ix-1,iy))+&
                     bvort(ix,iy)*vort(ix,Ny1+Ny2+1)
       else
         dvort(ix,iy)=vort(ix,iy)-(Cx*u(ix,iy)/4.d0)*&
                     (vort(ix+1,iy)-vort(ix-1, iy))+&
                     (Fox*Pr/2.d0)*(vort(ix+1,iy)+vort(ix-1, iy)-2.d0*vort(ix,iy))&
                    -(Cx*Ri/4.d0)*(temp(ix+1,iy)-temp(ix-1,iy))
      
       endif
     enddo
   endif
 enddo
!========Resolution de la vorticite vort(n+1)=======
 do ix=2,Nx1+Nx2
   if (ix.le.Nx1+1) then
      call TDMA(avort(ix,(Ny1+2):(Ny1+Ny2)), bvort(ix,(Ny1+2):(Ny1+Ny2)),&
               cvort(ix,(Ny1+2):(Ny1+Ny2)), dvort(ix,(Ny1+2):(Ny1+Ny2)),&
               vort(ix,(Ny1+2):(Ny1+Ny2)),Ny2-1)                      
   else
       call TDMA(avort(ix,2:(Ny1+Ny2)), bvort(ix,2:(Ny1+Ny2)),&
               cvort(ix,2:(Ny1+Ny2)), dvort(ix,2:(Ny1+Ny2)),&
               vort(ix,2:(Ny1+Ny2)),Ny1+Ny2-1)                      
   endif
 enddo
 
!==============RESOLUTION DE PSI PAR LA METHODE SoR===========
  call bc(temp, psi, vort, u, v)  ! Les conditions aux limites
 comp=0
 diff=1000.d0
 do while ((diff.ge.error).and.(comp.le.Nmax))
   copyPsi = psi
   do iy=2,Ny1+Ny2
     if (iy.le.Ny1+1) then
       do ix=2,Nx2  
          psi(ix+Nx1, iy)=(1-gama)*psi(ix+Nx1,iy) + gama/(2*(1+(dx/dy)**2))*&
           (psi(ix+Nx1+1, iy) + psi(ix+Nx1-1, iy) + ((dx/dy)**2)*psi(ix+Nx1, iy+1)+&
           ((dx/dy)**2)*psi(ix+Nx1, iy-1) - ((dx)**2)*vort(ix+Nx1,iy))                  
       enddo
     else
       do ix=2,Nx1+Nx2     
          psi(ix, iy) = (1-gama)*psi(ix, iy) + gama/(2*(1+(dx/dy)**2))*&
            (psi(ix+1, iy) + psi(ix-1, iy) + ((dx/dy)**2)*psi(ix, iy+1)+&
            ((dx/dy)**2)*psi(ix, iy-1) - ((dx)**2)*vort(ix, iy))                          
       enddo
     endif
   enddo
    diff=maxval(abs(psi-copyPsi))
    comp=comp+1  
 enddo
 write(*,*) 'time step ==>' , istep, 'Nb ite SoR :',  comp, 'Erreur SoR :',  diff
 if (comp.eq.Nmax) then
   write(*,*) 'Nombre maximal d iteration depasse!'
 endif

 !============ actualisation de vitesse =========

 do iy=2,Ny1+Ny2
   if (iy.le.Ny1+1) then
     do ix=2,Nx2     
      u(ix+Nx1,iy)=(psi(ix+Nx1,iy+1)-psi(ix+Nx1,iy-1))/(2.d0*dy)
      v(ix+Nx1,iy)=-(psi(ix+Nx1+1,iy)-psi(ix+Nx1-1,iy))/(2.d0*dx)                    
     enddo
   else
     do ix=2,Nx1+Nx2     
      u(ix,iy)=(psi(ix,iy+1)-psi(ix,iy-1))/(2.d0*dy)
      v(ix,iy)=-(psi(ix+1,iy)-psi(ix-1,iy))/(2.d0*dx) 
     enddo
   endif
  enddo
 !==============================================

 !Ecriture des resultats
 !---------------------
 if (mod(istep,fr).eq.0) then
  base='velo'
  write(int2char, 100) istep
  filename= trim(adjustl(base))//trim(adjustl(int2char))
  open(unit=istep+300, file =filename)
  do ix = 1, Nx1+Nx2+1
    do iy=1,Ny1+Ny2+1
      write(300+istep,*) xx(ix),yy(iy), u(ix, iy), v(ix, iy), psi(ix, iy), vort(ix, iy), temp(ix, iy)
    enddo
    write(300+istep,*)  '   '
  enddo
  close(istep+300)
  endif
100 FORMAT(' ',I5.5)
 !==== Ecriture des fichiers pour lancer une suite de calcul ==
 if (istep.eq.Nstep) then
  open(unit=Nstep+500, file ='res_suite')
  if(Nstep.gt.Nstep_old) then
     write(500+Nstep,*) Nstep
  else
     stop
  endif
  do ix = 1, Nx1+Nx2+1
     do iy=1,Ny1+Ny2+1
        write(500+Nstep,*) xx(ix),yy(iy), u(ix, iy), v(ix, iy), psi(ix, iy), vort(ix, iy),temp(ix, iy)
     enddo
     write(500+istep,*)  '   '
  enddo
  close(Nstep+500)
  endif
 !=== Calcul du coefficient de frottement sur les parois B2 B3 B4 =============
 if (istep.eq.Nstep) then
  open(unit=Nstep+8000, file ='coeff_frott.dat')
  write(8000+Nstep,*) '#1)- X     2)- Cf'
  do ix = 1, Nx1+Nx2+1
     if(ix.ge.Nx1+1) then
         write(8000+Nstep,*) xx(ix) , 2*u(ix,2)/(dy*Re)
     else
         write(8000+Nstep,*) xx(ix) , 2*u(ix,Ny1+2)/(dy*Re) 
     endif
   enddo
  close(Nstep+8000)
  endif
 !========profils de vitesse a x=1, 2, 3, 4, 5, 6, 7, 8 ====================
  !=====Profil a x=H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8003, file ='profil1H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-1.d0).le.1.d-5) then
        write(8003+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8003)
  endif
  !=====Profil a x=2H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8004, file ='profil2H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-2.d0).le.1.d-5) then
        write(8004+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8004)
  endif
  !=====Profil a x=3H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8005, file ='profil3H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-3.d0).le.1.d-5) then
        write(8005+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8005)
  endif
  !=====Profil a x=4H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8006, file ='profil4H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-4.d0).le.1.d-5) then
        write(8006+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8006)
  endif
  !=====Profil a x=5H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8007, file ='profil5H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-5.d0).le.1.d-5) then
        write(8007+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8007)
  endif
  !=====Profil a x=6H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8008, file ='profil6H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-6.d0).le.1.d-5) then
        write(8008+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8008)
  endif
   !=====Profil a x=7H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8009, file ='profil7H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-7.d0).le.1.d-5) then
        write(8009+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8009)
  endif
    !=====Profil a x=8H
 if (istep.eq.Nstep) then
  open(unit=Nstep+8010, file ='profil8H.dat')
  do iy=1,Ny1+Ny2+1 
     do ix=1,Nx2+1
       if(abs(xx(ix+Nx1)-xx(1+Nx1)-8.d0).le.1.d-5) then
        write(8010+Nstep,*) yy(iy), u(ix+Nx1, iy), v(ix+Nx1, iy), psi(ix+Nx1, iy), vort(ix+Nx1, iy), temp(ix+Nx1, iy)
       endif
     enddo
  enddo
  close(Nstep+8010)
  endif
 
 !====================================================================
 !=== Calcul du Nombre de Nusselt sur les parois B2 B3 B4=============
 if (istep.eq.Nstep) then
  open(unit=Nstep+8001, file ='nusselt_profileX.dat')
  open(unit=Nstep+8002, file ='nusselt_profileY.dat')
  write(8001+Nstep,*) '#1)- X     2)- Nu'
  write(8002+Nstep,*) '#1)- X     2)- Nu'
  do ix = 1, Nx1+Nx2+1
     if(ix.gt.Nx1+1) then
         write(8001+Nstep,*) xx(ix) , (temp(ix,2)-temp(ix,1))/dy
     else
         write(8001+Nstep,*) xx(ix) , (temp(ix,Ny1+2)-temp(ix,Ny1+1))/dy 
     endif
   enddo
  do iy = 1, Ny1+1
      write(8002+Nstep,*) yy(iy) , (temp(Nx1+2,iy)-temp(Nx1+1,iy))/dx
  enddo
  close(Nstep+8001)
  close(Nstep+8002)
  endif
 !=====================================================
 !====impression des variables au niveau de la sonde ===
 write(11,*) istep , u(ixs,iys), v(ixs,iys), psi(ixs,iys), vort(ixs,iys), temp(ixs,iys)
 !================================================
enddo ! fin de la boucle sur le temps
endif
!=================================================
!*********************************************
call system('gnuplot script.GNU')
call system('gnuplot profile.GNU')
open(unit=8, file = 'mesh') 
do ix = 1, Nx1+Nx2+1
  if(ix.ge.Nx1+1) then
    do iy=1,Ny1+Ny2+1
       write(8,*) xx(ix) , yy(iy)
    enddo
  else
    do iy=1, Ny2+1
       write(8,*) xx(ix) , yy(iy+Ny1)
    enddo
 endif
enddo
close(8)
close(11)

11 format(20(1E14.7, ' '))

write(*,*) 'FourierX', Fox
write(*,*) 'FourierY', Foy
write(*,*) 'CourantX', Cx
write(*,*) 'CourantY', Cy
end 


