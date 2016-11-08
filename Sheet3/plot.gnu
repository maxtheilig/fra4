f(x) = 0.0185695;
set xrange [100:10000]
set yrange [0.0:0.04]
set terminal png
set output 'steps_result.png'
plot 'steps_result.dat' using 1:2, f(x)
exit
