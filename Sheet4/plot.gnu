set terminal png
set output 'Wfluc.png'
f(x) = 1
set title "Fluctuation of W"
set xlabel "lambda_1"
set ylabel "<W^2>/<W>^2"
unset key
set xrange [-11:2]
set yrange [0:2]
plot 'Wfluc.txt' using 1:2 with points pointtype 5, f(x)
exit
