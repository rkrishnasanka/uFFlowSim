#!/bin/sh
# X=0.0 # Reference X coordinate in top left corner
# Y=0.0 # Reference Y coordinate in top left corner
# W=1.0 # Centered
# H=0.1 # Y delta from top; doesn't work is set smaller
# D=0.1 # Quantity of fluid coming from source
# U=0.0 # X direction of flow magnitude
# V=3.0 # Y direction of flow magnitude

cd ..
./curved-boundaries components/figure5.csv
rm -f output5.gif
ffmpeg -loglevel panic -i Frame00%03d.png output5.gif
rm *.png
echo "Done."
