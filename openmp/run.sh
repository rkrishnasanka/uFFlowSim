#!/bin/sh

# ./curved-boundaries figure1.csv
# rm output1.gif
# ffmpeg -loglevel panic -i Frame00%03d.png output1.gif
# rm *.png
# echo "Done."

for f in ./components/*.csv;

do
    cores=1

    while [ $cores -le 8 ]
    do
        export OMP_NUM_THREADS=$cores
        echo "Running File: $f, Cores: $OMP_NUM_THREADS" 
        ./curved-boundaries $f
        ffmpeg -loglevel panic -i Frame00%03d.png output.gif
        rm *.png
        mv output.gif 'results/output_'$(basename $f)'_'$cores'.gif';
        echo "=========Done=========="
        ((cores*=2))
    done


done