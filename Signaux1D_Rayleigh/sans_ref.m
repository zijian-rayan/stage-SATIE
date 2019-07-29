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
imageoy = real(ifft2(ifftshift(imaxx)));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% direction z
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
imageoz = real(ifft2(ifftshift(imaxx)));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% direction x
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
imageox = real(ifft2(ifftshift(imaxx)));





%%
figure(3)
subplot(221)
    plot(t,imageox(200,:))
    xlabel('Time (s)')
    ylabel('Vibration Vx')
    title('Vx(t) for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])
subplot(222)
    plot(t,imageoy(200,:))
    xlabel('Time (s)')
    ylabel('Vibration Vy')
    title('Vy(t) for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])
subplot(223)
    plot(t,imageoz(200,:))
    xlabel('Time (s)')
    ylabel('Vibration Vz')
    title('Vz(t) for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])

%All
subplot(224)
    plot(t,imageox(200,:),t,imageoy(200,:),t,imageoz(200,:))
    xlabel('Time (s)')
    ylabel('Vibrations Vx, Vy, Vz')
    title('V for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])

%Ellipticity:
%theoretical ellipso H/V    
qt=(2765/3388)^2;
hvt=2*(sqrt(1-qt)./(2-qt)) %.866


figure(4)

int=120:190;
subplot(211)
    plot(t(int),imageox(200,int),t(int),imageoz(200,int))
    xlabel('Time (s)')
    ylabel('Vibrations Vx, Vz')
    title('V for x=44.6mm')
    %axis([0 5e-5 -5e-3 5e-3])
subplot(212)
mmax=max([imageox(200,int),imageoz(200,int)]');
mvvx=imageox(200,int)-mean(imageox(200,int));
mvvz=imageoz(200,int)-mean(imageoz(200,int));
mmax=max([mvvx,mvvz]');

    [zet]=elliptic_plot(mmax*hvt,mmax,0,0,0);
    [ze]=elliptic_plot(mmax*.5,mmax,0,0,0);
    plot(imageox(200,int),imageoz(200,int))
    xlabel('Vx')
    ylabel('Vz')
    title('Ellipsometry for x=44.6mm')
    axis([-5e-3 5e-3 -5e-3 5e-3])
    axis equal
%theoretical ellipso H/V    
q=(2765/3388)^2;
hv=2*(sqrt(1-q)./(2-q))
%%


