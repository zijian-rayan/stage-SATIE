
%%
%
%
%
%       Cette section pr¨¦sente une comparaison de bout en bout 
%d      u filtrage du domaine fr¨¦quentiel dans trois directions.
%
%
%%
clear all;
close all;
clc;
importfiles=0;
load Rayleigh3D.mat
%%
%Caution: z on propa axis, x normal to surf.... for ALU
vvx=flipud(vz);
vvy=flipud(vy);
vvz=flipud(vx);
[nx nt]=size(vx);
%goto complex?
% cv=squeeze(vvy(:,5,1:600))+sqrt(-1)*squeeze(vvz(:,5,1:600));
%Caution: z on propa axis, x normal to surf.... for ALU
cv=sqrt(-1)*squeeze(vvx(:,1:600))+squeeze(vvz(:,1:600));

dx=abs(min(diff(XYZ(:,3))));
x=(0:nx-1)*dx;
dt=min(diff(t))
f=((0:nt-1)/nt-.5)/dt;

%Compare 3 directions on all FOV:

mmax=max([vvx(:);vvy(:);vvz(:)])/10;

  load config_file
  load signal_600_1.mat
       
%%  direction y
sigma = 1.6;
RR =fftshift(abs(fft2(vvy)));% fft
RR = wiener2(RR,[5 10]); %filtre wiener2D
thresh=[0.10,0.11];
% thresh=[0.01,0.17];
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
se = strel('disk',5);

f=imdilate(abs(f),se);
f=imerode(f,se);
    f(:,1180:1320)=0;
    f(1:50,:)=0;
    f(387:437,:)=0;
imaxx = fftshift((fft2(vvy)));
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
gausFilter = fspecial('gaussian', [5,5], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(vvy)));
imaxx=imaxx.*gaus;
imageo2 = real(ifft2(ifftshift(imaxx)));
imaxor = imageo2-imageo;
figure(1)
subplot(323)
    imagesc(t,x,vvy,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('y vibration')
subplot(324)
    imagesc(t,x,imageo2,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('y vibration')
mean_y =mean(max(imageo2))
figure(2)
subplot(132)
mesh(imageo2)
title('mesh(y)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% direction z
figure(1)
sigma = 1.6;
RR =fftshift(abs(fft2(vvz)));% fft
RR = wiener2(RR,[5 10]); %filtre wiener2D
thresh=[0.11,0.37];
% thresh=[0.01,0.17];
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);

se = strel('disk',5);

f=imdilate(abs(f),se);
f=imerode(f,se);
    f(:,1180:1320)=0;
    f(1:50,:)=0;
    f(387:437,:)=0;
imaxx = fftshift((fft2(vvz)));
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
gausFilter = fspecial('gaussian', [5,5], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(vvz)));
imaxx=imaxx.*gaus;
imageo2 = real(ifft2(ifftshift(imaxx)));
imaxor = imageo2-imageo;
subplot(325)
    imagesc(t,x,vvz,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('z vibration')
subplot(326)
    imagesc(t,x,imageo2,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('z vibration')
mean_z = mean(max(imageo2))

figure(2)
subplot(133)
mesh(imageo2)
title('mesh(z)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% direction x
figure(1)
sigma = 1.6;
RR =fftshift(abs(fft2(vvx)));% fft
RR = wiener2(RR,[5 10]); %filtre wiener2D
thresh=[0.11,0.37];
% thresh=[0.01,0.17];
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);

 se = strel('disk',5);

f=imdilate(abs(f),se);
f=imerode(f,se);
    f(:,1180:1320)=0;
    f(1:50,:)=0;
    f(387:437,:)=0;
imaxx = fftshift((fft2(vvx)));
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
gausFilter = fspecial('gaussian', [5,5], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(vvx)));
imaxx=imaxx.*gaus;
imageo2 = real(ifft2(ifftshift(imaxx)));
imaxor = imageo2-imageo;
subplot(321)
    imagesc(t,x,vvx,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('x vibration')
subplot(322)
    imagesc(t,x,imageo2,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('x vibration')
mean_x = mean(max(imageo2))
figure(2)
subplot(131)
mesh(imageo2)
title('mesh(x)')


%%



