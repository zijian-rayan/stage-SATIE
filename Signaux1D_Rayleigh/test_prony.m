%Prony check noise
%%%One pole
clc
clear
close all
n=2000;
dx=.1e-3;
x=(0:n-1)*dx;
k=200+sqrt(-1)*20;
A=10;
s=A*exp(sqrt(-1)*k*x);
figure(1)
subplot(231)
plot(x,real(s));
title('Signal origine')
%%
%Prony
subplot(232)
nb=0;
na=1;
[bn,an]=prony(s,nb,na);
[rn1,pn1,kn1]=residuez(bn,an);
%p=exp(jkdx)
kn1=-sqrt(-1)*log(pn1)/dx;
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);%³å»÷ÏìÓ¦
plot(x,real(s),x,real(s),'k',xn1,real(sn1),'r')
title('Impz 0 1')
subplot(233)
%%%Same with noise
noise=rand(1,n);
noise=noise-mean(noise);
sn=s+noise;
plot(real(sn))
title('Ajouter du bruit')
%%
subplot(234)
nb=0;
na=1;
[bn,an]=prony(sn,nb,na);
[rn1,pn1,kn1]=residuez(bn,an);
%p=exp(jkdx)
kn1=-sqrt(-1)*log(pn1)/dx;
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);%³å»÷ÏìÓ¦
plot(x,real(s),'b',xn1,real(sn1),'r')
legend('s','sn','snl');
title('0 1 estimation')
%%
%Same more poles
subplot(235)
noise=rand(1,n);
sn=s+noise;
plot(real(sn))
title('Signal contenant du bruit')
%%
subplot(236)
nb=9;
na=10;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('9 10 estimation')
%%
%-->quadratic error
error=sum(abs(sall-sn')).^2;
cs_all_sn=(xcorr(sall,conj(sn')));
%loop
nn=20;
error =zeros(nn,1)
for i=1:nn;
    [bn,an]=prony(sn,i-1,i);
    [sall,xall]=impz(bn,an,n,1/dx);
    error(i)=sum(abs(sall-sn')).^2;
    cs_all_sn(i)=max(xcorr(sall,conj(sn')));

end
figure
subplot(331)
nb=0;
na=1;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('0 1 estimation')
%%
subplot(332)
nb=1;
na=10;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('1 10 estimation')
%%
subplot(333)
nb=1;
na=100;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('1 100 estimation')
%%
subplot(334)
nb=2;
na=1;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('2 1 estimation')
%%
subplot(335)
nb=10;
na=1;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('10 1 estimation')
%%
subplot(336)
nb=100;
na=1;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('100 1 estimation')
%%
subplot(337)
nb=9;
na=10;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('9 10 estimation')
%%
subplot(338)
nb=10;
na=20;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('10 20 estimation')
%%
subplot(339)
nb=300;
na=300;
[bn,an]=prony(sn,nb,na);
[rn,pn,kn]=residuez(bn,an);
%p=exp(jkdx)
kn=-sqrt(-1)*log(pn)/dx;
%--> residue rn good as fast and fair selection criteria
% compare ident with original
[sall,xall]=impz(bn,an,n,1/dx);
[vmax,imax]=fmax2(abs(rn));
[smain,xmain]=impz(rn(imax(1)),[1 -pn(imax(1))],n,1/dx);
[sn1,xn1]=impz(rn1,[1 -pn1],n,1/dx);
plot(x,real(s),'k',xall,real(sall),'r',xmain,real(smain),'c')
legend('s','sall','snmain');
title('300 300 estimation')
%%