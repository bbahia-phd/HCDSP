



function crg_mat_wigb(dx,dy,dz;
    dt=0.004,
    nt=100,n1=20,npanel=3,
    irx=10,iry=10,
    xtickpos = [16;42;67;93],
    xticklab = ([x[5];x[10];x[15];x[20]] .- 1) .* 10,
    ytickpos=1:50:99,
    FontSize=16,
    sc=0.5,
    amx=1.0,
    fname="noisy5d3c",
    w=10,
    h=15)
t = collect((0:nt-1).*dt);

tmpx = zeros(nt,n1,npanel);
tmpy = zeros(nt,n1,npanel);
tmpz = zeros(nt,n1,npanel);
j = 1;
for i in (10,15,20)
tmpx[:,:,j] .= dx[:,:,i,10,10];
tmpy[:,:,j] .= dy[:,:,i,10,10];
tmpz[:,:,j] .= dz[:,:,i,10,10];
j += 1;
end

tmpx[1,:,:] .= 0;
tmpy[1,:,:] .= 0;
tmpz[1,:,:] .= 0;
tt,ty,tx = size(tmpx);
x = zeros(tx*ty);

c0=6;
k=1; c=0;
for ix = 1:tx
c += c0;
for iy = 1:ty
x[k] = (ix-1)*(ty-1)+iy+c;
k += 1;
end
end

Dx = reshape(tmpx,nt,:);
Dy = reshape(tmpy,nt,:);
Dz = reshape(tmpz,nt,:);


mat"close all"
mat"h3 = tight_subplot(3,1,0.05,[0.1 0.07],[0.1 0.01])"
mat"subplot(311)"
mat"wigb($(Dx),$(sc),$x,$t,$(amx))"
mat"title({'{(a)} Ux'},'FontSize',$FontSize,'Interpreter','latex')"
mat"ylabel({'t(s)'},'FontSize',$FontSize,'Interpreter','latex')"
mat"xlabel({'sx (m)'},'FontSize',$FontSize,'Interpreter','latex')"
mat"set(gca,'TickLabelInterpreter','latex')"
mat"set(gca,'FontSize',$FontSize)"
mat"xticks($xtickpos)"
mat"xticklabels($xticklab)"
mat"yticks([0 0.2 0.396])"
mat"yticklabels([0 0.2 0.4])"

mat"subplot(312)"
mat"wigb($(Dy),$(sc),$x,$t,$(amx))"
mat"title({'{(b)} Uy'},'FontSize',$FontSize,'Interpreter','latex')"
mat"ylabel({'t(s)'},'FontSize',$FontSize,'Interpreter','latex')"
mat"xlabel({'sx (m)'},'FontSize',$FontSize,'Interpreter','latex')"
mat"set(gca,'TickLabelInterpreter','latex')"
mat"set(gca,'FontSize',$FontSize)"
mat"xticks($xtickpos)"
mat"xticklabels($xticklab)"
mat"yticks([0 0.2 0.396])"
mat"yticklabels([0 0.2 0.4])"

mat"subplot(313)"
mat"wigb($(Dz),$(sc),$x,$t,$(amx))"
mat"title({'{(c)} Uz'},'FontSize',$FontSize,'Interpreter','latex')"
mat"ylabel({'t (s)'},'FontSize',$FontSize,'Interpreter','latex')"
mat"xlabel({'sx (m)'},'FontSize',$FontSize,'Interpreter','latex')"
mat"set(gca,'TickLabelInterpreter','latex')"
mat"set(gca,'FontSize',$FontSize)"
mat"xticks($xtickpos)"
mat"xticklabels($xticklab)"
mat"yticks([0 0.2 0.396])"
mat"yticklabels([0 0.2 0.4])"

mat"set(gcf,'PaperOrientation','landscape')"
mat"set(gcf,'PaperSize',[$w $h])"
mat"print(gcf,$(fname),'-dpdf','-fillpage')"
#    mat"exportgraphics(gcf,$(fname),'ContentType','vector')"
end