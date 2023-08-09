step='step.geo'
set terminal pngcairo size  1024,768
############# U profiles with the step ##################
set output 'profile_velo_step.png'
set size ratio -1
set grid
set key off
set xtics 1
set title 'profils des vitesses après la marche'
p step w l lw 1.5 lc 0,\
'profil1H.dat' u ($2+10+1):1 notitle w lp pt 3 ps 0.3,\
'profil2H.dat' u ($2+10+2):1 notitle w lp pt 3 ps 0.3,\
'profil3H.dat' u ($2+10+3):1 notitle w lp pt 3 ps 0.3,\
'profil4H.dat' u ($2+10+4):1 notitle w lp pt 3 ps 0.3,\
'profil5H.dat' u ($2+10+5):1 notitle w lp pt 3 ps 0.3,\
'profil6H.dat' u ($2+10+6):1 notitle w lp pt 3 ps 0.3,\
'profil7H.dat' u ($2+10+7):1 notitle w lp pt 3 ps 0.3,\
'profil8H.dat' u ($2+10+8):1 notitle w lp pt 3 ps 0.3
###################V profiles ###########################
reset 
set terminal pngcairo size  1024,768
set output 'V.png'
set grid
set title 'Profils de la deuxième compotante de vitesse V'
set xlabel 'V/U0'
set ylabel 'Y/H'
set key t l 
p 'profil1H.dat' u 3:1 t 'x=1H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil2H.dat' u 3:1 t 'x=2H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil3H.dat' u 3:1 t 'x=3H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil4H.dat' u 3:1 t 'x=4H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil5H.dat' u 3:1 t 'x=5H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil6H.dat' u 3:1 t 'x=6H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil7H.dat' u 3:1 t 'x=7H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil8H.dat' u 3:1 t 'x=8H' w lp pt 3 ps 0.3 lw 1.5
##################### U profiles  #########################
reset 
set terminal pngcairo size  1024,768
set output 'U.png'
set grid
set title 'Profils de la première compotante de vitesse U'
set xlabel 'U/U0'
set ylabel 'Y/H'
set key t l 
p 'profil1H.dat' u 2:1 t 'x=1H' w lp pt 3 ps 0.3 lw 1.2,\
  'profil2H.dat' u 2:1 t 'x=2H' w lp pt 3 ps 0.3 lw 1.2,\
  'profil3H.dat' u 2:1 t 'x=3H' w lp pt 3 ps 0.3 lw 1.2,\
  'profil4H.dat' u 2:1 t 'x=4H' w lp pt 3 ps 0.3 lw 1.2,\
  'profil5H.dat' u 2:1 t 'x=5H' w lp pt 3 ps 0.3 lw 1.2,\
  'profil6H.dat' u 2:1 t 'x=6H' w lp pt 3 ps 0.3 lw 1.2,\
  'profil7H.dat' u 2:1 t 'x=7H' w lp pt 3 ps 0.3 lw 1.2,\
  'profil8H.dat' u 2:1 t 'x=8H' w lp pt 3 ps 0.3 lw 1.2
############### Zoom on U ##################################
reset 
set terminal pngcairo size  1024,768
set output 'U_zoom.png'
set grid
set title 'Profils de la première compotante de vitesse U'
set yrange[0:1.8]
set xlabel 'U/U0'
set ylabel 'Y/H'
set key t l 
p 'profil1H.dat' u 2:1 t 'x=1H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil2H.dat' u 2:1 t 'x=2H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil3H.dat' u 2:1 t 'x=3H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil4H.dat' u 2:1 t 'x=4H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil5H.dat' u 2:1 t 'x=5H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil6H.dat' u 2:1 t 'x=6H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil7H.dat' u 2:1 t 'x=7H' w lp pt 3 ps 0.3 lw 1.5,\
  'profil8H.dat' u 2:1 t 'x=8H' w lp pt 3 ps 0.3 lw 1.5
################ Skin friction coefficient ##########################
reset 
set terminal pngcairo size 1024,768  
set xtics 2
set size ratio 0.5
set output 'Cf.png'
set grid
set title 'Profil du coefficient de frottement Cf'
set xlabel 'X/H'
set ylabel 'Cf'
set key off
p step u 1:($2/520) w l lc 0 lw 2,\
'coeff_frott.dat' u 1:2  w p pt 3 ps 0.1 lw 2 lc 1
reset 

############# Nusselt profile  X ###########################
set terminal pngcairo size  1024,768
set output 'Nu.png'
set multiplot layout 2,1        # engage multiplot mode
set grid
set title 'Profil du nombre de Nusselt Nu suivant X'
set xlabel 'X/H'
set ylabel 'Nu'
set xtics 2
set key off
p step u 1:($2/4) w l lc 0 lw 2,\
'nusselt_profileX.dat' u 1:2 notitle w lp pt 3 ps 0.3 lw 2 lc 1 
############# Nusselt profile  Y ###########################
set grid
set xtics 1
set title 'Profil du nombre de Nusselt Nu suivant la hauteur de la marche'
set ylabel 'Y/H'
set xlabel 'Nu'
set key off
p 'nusselt_profileY.dat' u 2:1 notitle  w lp pt 3 ps 0.3 lw 2 lc 1 
unset multiplot

#########################################################
############################
file='res_suite'
###########################

reset 

#----------U -------------
set terminal pngcairo size 1024,768
set output 'vitesseU.png'
set title 'Profil de la première composante de  vitesse U'
set key off
set size ratio -1
set view map
set contour base
set cntrparam levels incremental -10,0.05,10
unset surface
set pm3d map
set pm3d interpolate 3,3
set surface
set cbrange [*:*]
set cblabel 'U/U0'
set palette defined
splot file  u 1:2:3 notitle w l , file u 1:2:3 notitle

reset 
#----------V -------------
set terminal pngcairo size 1024,768
set output 'vitesseV.png'
set title 'Profil de la deuxième composante de  vitesse V'
set pm3d map
set pm3d interpolate 3,3
set cbrange [*:*]
set cblabel 'V/U0'
set key off
set size ratio -1
set palette defined
splot  file  every ::1 u 1:2:4
reset 
#----------T -------------
set terminal pngcairo size 1024,768
set output 'temperature.png'
set title 'Profil de la temperature adimensionnée'
set pm3d map
set pm3d interpolate 3,3
set cbrange [0:1]
set cblabel '(T-T0)/(T1-T0)'
set key off
set size ratio -1
set palette defined
splot  file  every ::1 u 1:2:7

reset 
#----------U contour -------------
set terminal pngcairo size 1024,768
set output 'vitesse_contour.png'
set title 'iso-contour de la vitesse'
set key off
set size ratio -1
set view map
set contour base
set cntrparam levels incremental -10,0.05,10
unset surface
set cbrange [*:*]
set cblabel 'U/U0'
set palette defined
splot file  u 1:2:(sqrt($3**2+$4**2)) notitle w l lt -1 
reset 

#----------vector -------------
set terminal pngcairo size 1024,768
set output 'vector.png'
set size ratio 0.5
set title 'Champ de vitesse juste après la marche'
set cbrange [0:1]
set xrange[9.5:18]
set yrange[0:1.1]
set key off
set palette defined
plot  file  u 1:2:3:4 w vector lc 0 



system('mkdir Figures')
system('mv *png Figures')





