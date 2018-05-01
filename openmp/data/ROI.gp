set xlabel 'Candidate Functions' offset 0, 4
set ylabel 'Percentage of Runtime'
set title 'ROI Analysis'
set boxwidth 0.2
set style fill solid
total=95.666634
set xtics rotate
set xtics font "Times-New-Roman, 9"
set terminal png
set output 'ROI.png'
plot 'ROI.dat' u 1:($2/total):xtic(3) w boxes t ''
