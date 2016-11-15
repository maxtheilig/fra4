set terminal png
set output 'error.png'
set title "Expectation Value of x^2 with error bars"
set xlabel "lambda_1"
set ylabel "<x^2>"
unset key
set xrange [-11:2]
#set yrange [0:2]
set style line 2 lc rgb 'black' pt 7   # circle
plot 'testerror.txt' using 1:2:3 pt 7 ps 1 lc rgb 'black' with errorbars, '1/xsquared.txt' using 1:2 pt 7 ps 1
exit
