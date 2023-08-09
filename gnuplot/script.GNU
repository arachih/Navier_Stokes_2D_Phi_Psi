#set terminal gif
###################### Temperature
set terminal pngcairo size 800,800 
set pm3d map
set pm3d interpolate 3,3
#set cbrange [-0.5:1.2]
set cbrange [0:1]
set key off
set size ratio -1
set palette defined
filelist = system('ls velo*')
do for [file in filelist] {
    set output 'Temp.'.(file[2:*]).'.png'
#    plot  file  u 1:2:7  palette  ps 1 pt 5
#    plot  file  u 1:2:3:4  w vector lc 0
     splot  file  u 1:2:7  
#    p file u 1:2:7 w image notitle,for [i=1:505] -0.09<=x && x<=9.99 ? i*0.002-0.048 :1/0 notitle lc rgb "white" lw 4

}

#########Première composante de vitesse U
reset
set terminal pngcairo size 800,800 
set pm3d map
set pm3d interpolate 3,3
set cbrange [-0.5:1.2]
#set cbrange [0:1]
set key off
set size ratio -1
set palette defined
filelist = system('ls velo*')
do for [file in filelist] {
    set output 'vitX.'.(file[2:*]).'.png'
#    plot  file  u 1:2:7  palette  ps 1 pt 5
#    plot  file  u 1:2:3:4  w vector lc 0
     splot  file  u 1:2:3
#    p file u 1:2:7 w image notitle,for [i=1:505] -0.09<=x && x<=9.99 ? i*0.002-0.048 :1/0 notitle lc rgb "white" lw 4

}


system('convert -delay 5 -loop 0 Temp*.png temp.gif')
system('convert -delay 5 -loop 0 vitX*.png vitesseX.gif')
system('rm -rf Temp* vitX*')
