#! /bin/sh
#


TF='Times-Roman'
LF='Times-Italic'
HF='Helvetica'
HFb='Helvetica-Bold'

# Get labels
pslabel t="(a)" f=$HFb  size=15 nsub=0 > a.ps 
pslabel t="(b)" f=$HFb  size=15 nsub=0 > b.ps 
pslabel t="(c)" f=$HFb  size=15 nsub=0 > c.ps 

#######################################################################################
pdfcrop *.pdf
mv *-crop.pdf *.pdf

rm  *-crop.*
