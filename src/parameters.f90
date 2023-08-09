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
! T :  Temperature T(i,j)
! Vort:  Vorticity Vc(i,j)
! psi: Courant function psi(i,j)
! 
!====================================
MODULE Parametres

implicit none
!***************************************************************!
!	Variables Generales pour le Solveur Navier-Stokes			!
!***************************************************************!
!************************************
! PARAMETRES D'ENNTREE
!************************************
! Physique
double precision, parameter :: g=9.81d0        !constante de la gravite
double precision, parameter :: beta= 0.0034d0  ! coefficient de dilatation de l'air
double precision, parameter :: nu = 15.6d-6    ! viscosite cinématique de l'air 
double precision, parameter :: kappa =20.d-6   ! diffusite thermique de l'air 
! Conditions aux limites 
double precision, parameter :: Tc= 21.0d0+273.15d0 ! Temperature chaude
double precision, parameter :: Tf= 2.d0+273.15d0 ! Temperature froide
double precision, parameter :: U0=(2.d0/3.d0)*1.5*2*7.8d-3 ! Vitesse d'entree

!===============
! parametres geometrique  
!   ___________________________________
!   |                                  |
!   |Ly2                               |
!   |                                  |
!   |_________                         |
!             | Ly1                    | 
!       Lx1   |________________________|
!                    Lx2
!       
double precision, parameter :: Ly1=1.d0     ! Hauteur de la marche 
double precision, parameter :: Ly2=8.d0     ! Hauteur de l'Entree
double precision, parameter :: Lx1=10.d0    ! distance entre l'entree et nez de la marche
double precision, parameter :: Lx2=20.d0    ! distance entre le nez de la marche et la sortie

!Maillage
integer, parameter :: Nx1=100  ! nombre de noeuds entre l'entree et le nez de la marche
integer, parameter :: Ny1=10   ! nombre de noeuds au niveau de la marche (H)

! parametre de calcul
double precision, parameter :: dtref= 8.d-1   ! Pas de temps
double precision, parameter :: gama = 1.725d0   ! coefficient de relaxation
integer, parameter :: Nstep  = 100  !nombre de pas de temps de la simulation
integer, parameter :: Nmax = 10000 !nombre maximal pour tester la convergence en SoR
integer, parameter :: fr = 100       !frequency of output result
integer, parameter :: isuite =0       ! 0: pas de suite, 1: Avec une suite
double precision, parameter :: error = 5.d-6

! Nombres adimensionnes physique
!--------------------------------
double precision, parameter ::Pr= 0.78d0  ! Nombre de prandtl
double precision, parameter ::Re= 50.d0 ! Nombre de Reynolds
double precision, parameter ::Ri= 1.d0  ! Nombre de Richardson

!sonde de controle
!------------------
double precision, parameter  :: Xsonde =12.d0   !Première coordonnée de la sonde
double precision, parameter  :: Ysonde =1.d0   !Deuxième coordonnée de la sonde

!=================================================
!   PARAMETRE A NE PAS CHANGER !
!=================================================

!== Adimentionnement par la heuteur de la marche

double precision, parameter :: Xout=Lx2/Ly1       ! la distance entre le nez de la marche et la sortie
double precision, parameter :: Yin = Ly2/Ly1      ! diam d'entrée
double precision, parameter :: Xstep= Lx1/Ly1     ! Xtep adim 
double precision, parameter :: Ystep= Ly1/Ly1     ! Hauteur de la marche
!double precision, parameter :: dt=dtref*U0/Ly1

!Pas de maillage
!----------------
double precision, parameter ::  dx = Xstep/Nx1
double precision, parameter ::  dy = Ystep/Ny1
integer , parameter :: Nx2=ceiling(Xout/dx)
integer , parameter :: Ny2=ceiling(Yin/dy)

! Nombres adimensionnes physique
!--------------------------------
!double precision, parameter ::Pr= nu/kappa  ! Nombre de prandtl
!double precision, parameter ::Re= U0*Ly1/nu ! Nombre de Reynolds
!double precision, parameter ::Ri= g*beta*(Tc-Tf)*Ly1/U0**2  ! Nombre de Richardson


double precision, parameter :: dt=dtref*Re*nu/Ly1**2
!Nombres de Fourier et Courant
!-------------------------------
double precision , parameter :: Fox = (1.d0/(Re*Pr))*dt/dx**2 ! Nombre de Fourier selon x 
double precision , parameter :: Foy = (1.d0/(Re*Pr))*dt/dy**2 ! Nombre de Fourier selon y
double precision , parameter :: Cx=dt/dx ! Nombre de Courant selon x 
double precision , parameter :: Cy=dt/dy ! Nombre de Courant selon y

END MODULE Parametres


