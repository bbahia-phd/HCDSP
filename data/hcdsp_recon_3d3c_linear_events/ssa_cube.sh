#! /bin/sh
#

TF='Times-Roman'
LF='Times-Italic'
HF='Helvetica'
HFb='Helvetica-Bold'

d1=0.1
d2=100
d3=100

bc=1
wc=-1

# Get labels
$CWPROOT/bin/pslabel t="(a)" f=$HFb  size=20 nsub=0 > a.ps 
$CWPROOT/bin/pslabel t="(b)" f=$HFb  size=20 nsub=0 > b.ps 
$CWPROOT/bin/pslabel t="(c)" f=$HFb  size=20 nsub=0 > c.ps 
$CWPROOT/bin/pslabel t="(d)" f=$HFb  size=20 nsub=0 > d.ps 
$CWPROOT/bin/pslabel t="(e)" f=$HFb  size=20 nsub=0 > e.ps 
$CWPROOT/bin/pslabel t="(f)" f=$HFb  size=20 nsub=0 > f.ps
$CWPROOT/bin/pslabel t="(g)" f=$HFb  size=20 nsub=0 > g.ps
$CWPROOT/bin/pslabel t="(h)" f=$HFb  size=20 nsub=0 > h.ps
$CWPROOT/bin/pslabel t="(i)" f=$HFb  size=20 nsub=0 > i.ps

#######################################################################################
# Generate ps files
$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_ssa_zx.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1 bclip=$bc wclip=$wc> tmp1.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_ssa_zy.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1 bclip=$bc wclip=$wc> tmp2.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_ssa_zz.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1 bclip=$bc wclip=$wc> tmp3.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_qssa_zx.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1  bclip=$bc wclip=$wc> tmp4.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_qssa_zy.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1 bclip=$bc wclip=$wc> tmp5.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_qssa_zz.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1 bclip=$bc wclip=$wc> tmp6.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_aqssa_zx.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1  bclip=$bc wclip=$wc> tmp7.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_aqssa_zy.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1 bclip=$bc wclip=$wc> tmp8.ps

$CWPROOT/bin/pscube < ./hcdsp_recon_3d3c_linear_events_aqssa_zz.bin n1=100 n2=40 n3=40 d1=0.004 f1=0 d2=10 f2=1 d3=10 f3=1 perc=99 size1=2.5 size2=1.5 size3=1.0 d1num=$d1 d2num=$d2 d3num=$d3 f2num=0 f3num=0 labelsize=15 label1='t(s)' label2='Offset x(m)' label3='Offset y(m)' interp=1 bclip=$bc wclip=$wc> tmp9.ps

# Edit Path to SU psmerge 
$CWPROOT/bin/psmerge in=tmp1.ps translate=-18.5,-19.5 in=a.ps translate=-17.5,-14.5 in=tmp2.ps translate=-15,-19.5 in=b.ps translate=-14,-14.5 in=tmp3.ps translate=-11.5,-19.5 in=c.ps translate=-10.5,-14.5 in=tmp4.ps translate=-18.5,-23.5 in=d.ps translate=-17.5,-18.5 in=tmp5.ps translate=-15,-23.5 in=e.ps translate=-14,-18.5 in=tmp6.ps translate=-11.5,-23.5 in=f.ps translate=-10.5,-18.5 in=tmp7.ps translate=-18.5,-27.5 in=g.ps translate=-17.5,-22.5 in=tmp8.ps translate=-15,-27.5 in=h.ps translate=-14,-22.5 in=tmp9.ps translate=-11.5,-27.5 in=i.ps translate=-10.5,-22.5>tmp.ps

ps2pdf -dPDFSETTINGS=/prepress -dEPSCrop tmp.ps

cp tmp.pdf PSCube.pdf 

pdfcrop PSCube.pdf

mv PSCube-crop.pdf PSCube.pdf

rm  tmp*
