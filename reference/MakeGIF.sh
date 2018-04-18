#!/bin/sh
rm ./output.gif
ffmpeg -i Frame00%03d.png output.gif
rm ./*.png
