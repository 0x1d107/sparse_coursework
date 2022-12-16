set term tikz;
set output 'comp.tex';
set xlabel 'N';
set ylabel 'время, с';
plot 'comp.data' using 1:2 title "LU-разложение" with lines, 'comp.data' using 1:3 title "BiCGStab" with lines;
