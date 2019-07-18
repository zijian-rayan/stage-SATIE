%% prony

clear all;
%close all;
clc;
load Rayleigh3D.mat
[nx nt]=size(vx);
dx=abs(min(diff(XYZ(:,3))));
x=(0:nx-1)*dx;
vvy=flipud(vy);
vvx=flipud(vz);
vvz=flipud(vx);
mmax=max([vvy(:)])/10;


for i=1:437
    ave(i)=mean(abs(vvy(i,:)));
end


filt_z=vvz;
filt_y1=vvy;
filt_x=vvx;


for i=1:437
        B = sort(unique(vvy(i,:)), 'descend');
        M = mean(vvy(i,:));
        for j=1:2500
            if abs(filt_y1(i,j))<B(20)
                filt_y1(i,j)=0;
            else
                if filt_y1(i,j)>0
                    filt_y1(i,j)=filt_y1(i,j)-B(20);
                else
                    filt_y1(i,j)=filt_y1(i,j)+B(20);
                end
            end
        end
end

for i=1:437
        B = sort(unique(vvz(i,:)), 'descend');
        M = mean(vvz(i,:));
        for j=1:2500
            filt_z(i,j)=filt_z(i,j)-M;
            if abs(filt_z(i,j))<B(20)
                filt_z(i,j)=0;
            else
                if filt_z(i,j)>0
                    filt_z(i,j)=filt_z(i,j)-B(20);
                else
                    filt_z(i,j)=filt_z(i,j)+B(20);
                end
            end
        end
end
for i=1:437
        B = sort(unique(vvx(i,:)), 'descend');
        M = mean(vvx(i,:));
        for j=1:2500
            filt_x(i,j)=filt_x(i,j)-M;
            if abs(filt_x(i,j))<B(20)
                filt_x(i,j)=0;
            else
                if filt_x(i,j)>0
                    filt_x(i,j)=filt_x(i,j)-B(20);
                else
                    filt_x(i,j)=filt_x(i,j)+B(20);
                end
            end
        end
end
%%
%figure(1)
sigma = 1.6;
RR =fftshift(abs(fft2(filt_y1)));% fft
RR = wiener2(RR,[5 10]); %filtre wiener2D
%thresh=[0.11,0.121];%y
thresh=[0.01,0.021];%new y
% thresh=[0.01,0.17];%xz
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
se = strel('disk',3);
se2 = strel('disk',5);
f=imdilate(abs(f),se);
f=imerode(f,se2);
    f(:,1180:1320)=0;
    f(1:50,:)=0;
    f(387:437,:)=0;
imaxx = fftshift((fft2(filt_y1)));
imaxx(f==0)=0;
imageo = real(ifft2(ifftshift(imaxx)));

bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
sigma = 1;
    fill(:,1180:1320)=0;
    fill(1:50,:)=0;
    fill(387:437,:)=0;
gausFilter = fspecial('gaussian', [3,3], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(filt_y1)));
imaxx=imaxx.*gaus;
imageoy = real(ifft2(ifftshift(imaxx)));

%%

sigma = 1.6;
RR =fftshift(abs(fft2(filt_z)));% fft
RR = wiener2(RR,[5 10]); %filtre wiener2D
%thresh=[0.11,0.121];%y
thresh=[0.11,0.37];%new 
% thresh=[0.01,0.17];%xz
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
se = strel('disk',3);
se2 = strel('disk',5);
f=imdilate(abs(f),se);
f=imerode(f,se2);
    f(:,1180:1320)=0;
    f(1:50,:)=0;
    f(387:437,:)=0;
imaxx = fftshift((fft2(filt_z)));
imaxx(f==0)=0;
imageo = real(ifft2(ifftshift(imaxx)));
bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
sigma = 1;
    fill(:,1180:1320)=0;
    fill(1:50,:)=0;
    fill(387:437,:)=0;
gausFilter = fspecial('gaussian', [3,3], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(filt_z)));
imaxx=imaxx.*gaus;
imageoz = real(ifft2(ifftshift(imaxx)));

%%

sigma = 1.6;
RR =fftshift(abs(fft2(filt_x)));% fft
RR = wiener2(RR,[5 10]); %filtre wiener2D
%thresh=[0.11,0.121];%y
thresh=[0.11,0.37];%new 
% thresh=[0.01,0.17];%xz
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
se = strel('disk',3);
se2 = strel('disk',5);
f=imdilate(abs(f),se);
f=imerode(f,se2);
    f(:,1180:1320)=0;
    f(1:50,:)=0;
    f(387:437,:)=0;
imaxx = fftshift((fft2(filt_x)));
imaxx(f==0)=0;
imageo = real(ifft2(ifftshift(imaxx)));

bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
sigma = 1;
    fill(:,1180:1320)=0;
    fill(1:50,:)=0;
    fill(387:437,:)=0;
gausFilter = fspecial('gaussian', [3,3], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(filt_x)));
imaxx=imaxx.*gaus;
imageox = real(ifft2(ifftshift(imaxx)));

%%
%$$$$$ prony  $$$$

%% Tentative de Prony 
Fs=10e6
dt=1/Fs
dx=0.1e-3
Fx=1/dx
%%
%s = real(imageox);
s = real(vvx);
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
figure(1)
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
%s = real(imageoy);
s = real(vvy);
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
   [r,p,k]=residue(b,a);%留数、极点以及直接项
	rr(i,1:length(r))=r';
	pp(i,1:length(p))=p';
   kk(i,1:length(k))=k';
end
figure(2)
subplot(221)
plot(fff,log(abs(pp(:,1)))*Fx,'.')
xlabel('Frequency (Hz)')
ylabel('Ksec')
title('Logarithme du p?le')
subplot(222)
%plot(ff,angle(pp(:,1))*Fx,'.')
imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold
subplot(223)
plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-')
%imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold,plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-'),hold
xlabel('Frequency (Hz)')
ylabel('K'' (rad/m)')
title('angle')
subplot(224)
%plot(ff,angle(pp(:,1))*Fx,'.')
VPest=abs(2*pi*fff'./fftshift(angle(pp(:,1))*Fx));
plot(fff,VPest,'r.-')
xlabel('Frequency (Hz)')
ylabel('Phase velocity (m/s)')
title('vitesse de phase')
%axis([-5e6 5e6 0 3000])
%%
%s = real(imageoz);
s = real(vvz);
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
s = real(imageox);
%s = real(vvx);
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
figure(4)
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
s = real(imageoy);
%s = real(vvy);
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
   [r,p,k]=residue(b,a);%留数、极点以及直接项
	rr(i,1:length(r))=r';
	pp(i,1:length(p))=p';
   kk(i,1:length(k))=k';
end
figure(5)
subplot(221)
plot(fff,log(abs(pp(:,1)))*Fx,'.')
xlabel('Frequency (Hz)')
ylabel('Ksec')
title('Logarithme du p?le')
subplot(222)
%plot(ff,angle(pp(:,1))*Fx,'.')
imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold
subplot(223)
plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-')
%imagesc(fff,kkk*2*pi,fftshift(abs(fft2(s)))),hold,plot(fftshift(fff),angle(pp(:,1))*Fx,'r.-'),hold
xlabel('Frequency (Hz)')
ylabel('K'' (rad/m)')
title('angle')
subplot(224)
%plot(ff,angle(pp(:,1))*Fx,'.')
VPest=abs(2*pi*fff'./fftshift(angle(pp(:,1))*Fx));
plot(fff,VPest,'r.-')
xlabel('Frequency (Hz)')
ylabel('Phase velocity (m/s)')
title('vitesse de phase')
%axis([-5e6 5e6 0 3000])
%%
s = real(imageoz);
%s = real(vvz);
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


