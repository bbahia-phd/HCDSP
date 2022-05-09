#!/bin/bash

for i in ./*.pdf; do
  pdfcrop "$i"
  mv *-crop.pdf $i
done
exit
