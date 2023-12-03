set terminal jpeg enhanced font 'arial,10' size 800,600
set output 'grafica_rdf.jpg'
set title 'Radial Distribution function g(r)'
set xlabel 'r[A]'
set ylabel 'g(r)'
plot 'rdf_0.2_temp_0.723.dat' using 1:2 
pause -1
set output 'grafica_coord.jpg'
