set terminal png
set output "autocorr.png"
set style data linespoints
plot "metro.txt" title "metropolis", "heat.txt" title "heatbath", "heat_or.txt" title "heatbath with overrelaxation"
