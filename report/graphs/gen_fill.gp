set term tikz;
set output 'fill.tex';
set xtics;
set xlabel "Заполненность матрицы";
set ytics;
set ylabel "Время работы, с";
set grid;
plot 'fill.data' using 1:2 title "LU-разложение" with lines, 'fill.data' using 1:3 title "BiCGStab" with lines;
