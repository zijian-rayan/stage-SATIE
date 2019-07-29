%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate 2D signal constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=10e6
dt=1/Fs
nt=2000
dx=0.1e-3
Fx=1/dx
nx=100
xx=(0:nx-1)*dx;
tt=(0:nt-1)*dt;
[t,x]=meshgrid(tt,xx);
%sin non attenuated
f0=1e6
k0=1e3
figure(1)
mesh(sin(2*pi*f0*t+2*pi*k0*x))
%%
%sin attenuated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Surface wave: 
%   *Wide band attenuated, non dispersive:
%   Start from space/freq half space S(x,omega):
%   => OK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
S=zeros(nx,nt/2);
f=(0:nt/2-1)/(nt/2)*Fs/2;
CR=2000;
Kp=2*pi*f/CR;
Ks=100
K=ones(nx,1)*(-Kp+sqrt(-1)*Ks);
Sp=exp(sqrt(-1)*K.*x(:,1:nt/2));
%s=exp(sqrt(-1)*K.*x(:,1:nt/2));
S=[Sp fliplr(conj(Sp))];
%mesh(real(ifft(S,[],4)))
%calc space-time signal:
s=real(ifft(S,[],2));
%%
mesh(s)
title('Signal original bas¨¦ sur un signal sinuso?dal')
%%
SS=fft(s,[],2);
freq=((0:nt-1)/(nt)-.5)*2*pi;
fff=((0:nt-1)/(nt)-.5)*Fs;
kkk=fliplr((0:nx-1)/(nx)-.5)*Fx;

%Prony identification:
% degree of den num important:
% 1 wave needs 2 conjugates poles
clear rr pp kk aa bb
for i=1:nt
%[nnum,dden]=prony(1,2,ss(:,i));
	[b,a]=prony(SS(:,i),1,3);
	aa(i,1:length(a))=a;
	bb(i,1:length(b))=b;
   [r,p,k]=residue(b,a);
	rr(i,1:length(r))=r';
	pp(i,1:length(p))=p';
   kk(i,1:length(k))=k';
end
figure(3)
subplot(221)
plot(fff,log(abs(pp(:,1)))*Fx,'.')
xlabel('Frequency (Hz)')
ylabel('Ksec')
subplot(222)
%plot(ff,angle(pp(:,1))*Fx,'.')
imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold
subplot(223)
plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-')
%imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold,plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-'),hold
xlabel('Frequency (Hz)')
ylabel('K'' (rad/m)')
subplot(224)
%plot(ff,angle(pp(:,1))*Fx,'.')
VPest=abs(2*pi*fff'./fftshift(angle(pp(:,1))*Fx));
plot(fff,VPest,'r.-')
xlabel('Frequency (Hz)')
ylabel('Phase velocity (m/s)')
%axis([-5e6 5e6 0 3000])
%%
figure(4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Surface wave: 
%   *Wide band attenuated, dispersive:
%   Start from space/freq half space S(x,omega):
%   => OK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=zeros(nx,nt/2);
f=(0:nt/2-1)/(nt/2)*Fs/2;
%CR=2000;
CR=4000*(atan(1-f/3e6)+pi/2)/pi;
Kp=2*pi*f./CR;
Ks=100
K=ones(nx,1)*(-Kp+sqrt(-1)*Ks);
%Sp=exp(sqrt(-1)*-Kp+sqrt(-1)*100.*x(:,1:nt/2))
Sp=exp(sqrt(-1)*K.*x(:,1:nt/2));
S=[Sp fliplr(conj(Sp))];
%mesh(abs(ifft(S)));
%calc space-time signal:
s=real(ifft(S,[],2));
%%
mesh(s)
title('Signal original bas¨¦ sur un signal tangent')
SS=fft(s,[],2);
freq=((0:nt-1)/(nt)-.5)*2*pi;
fff=((0:nt-1)/(nt)-.5)*Fs;
kkk=fliplr((0:nx-1)/(nx)-.5)*Fx;

%Prony identification:
% degree of den num important:
% 1 wave needs 2 conjugates poles
clear rr pp kk aa bb
for i=1:nt
%[nnum,dden]=prony(1,2,ss(:,i));
	[b,a]=prony(SS(:,i),1,2);
	aa(i,1:length(a))=a;
	bb(i,1:length(b))=b;
   [r,p,k]=residue(b,a);
	rr(i,1:length(r))=r';
	pp(i,1:length(p))=p';
   kk(i,1:length(k))=k';
end
figure(6)
subplot(221)
plot(fff,log(abs(pp(:,1)))*Fx,'.')
xlabel('Frequency (Hz)')
ylabel('Ksec')
subplot(222)
%plot(ff,angle(pp(:,1))*Fx,'.')
imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold
subplot(223)
plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-')
%imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold,plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-'),hold
xlabel('Frequency (Hz)')
ylabel('K'' (rad/m)')
subplot(224)
%plot(ff,angle(pp(:,1))*Fx,'.')
VPest=abs(2*pi*fff'./fftshift(angle(pp(:,1))*Fx));
plot(fff,VPest,'r.-')
xlabel('Frequency (Hz)')
ylabel('Phase velocity (m/s)')
%axis([-5e6 5e6 0 3000])

