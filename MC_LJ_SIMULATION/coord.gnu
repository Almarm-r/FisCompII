set view 60,30,1.0,1.5
set title 'coord'
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
splot 'coordinates.dat' using 1:2:3
pause -1
