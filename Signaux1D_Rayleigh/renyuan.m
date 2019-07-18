
% Checking Polytec import Matlab functions for svd files:
%(1)Prior use, install Polytec libraries:
%   - Polytec_Common_Runtime_171_64Bit_(1)
%   - Polytec_FileAccess_471_(2)


clear all;
close all;
clc;
importfiles=0;
if importfiles==1,
    
% opengl('save','software');
global pathcurrent;
pathcurrent = pwd;
%pathName = 'C:\Users\BILY_HOLLOW\Desktop\Stage_2015\DATA\PSV_3D\Marbre\Format_svd\';
pathName = 'D:\ÊµÑéÊÒÊµÏ°\Signaux1D_Rayleigh\';


%fileName = 'Scan_time_Marbre_g1_Alt_500KHZ.svd';
%fileName = 'Scan_time_Marbre_g2_Alt_500KHZ.svd';
%fileName = 'Scan_time_MarbreI1_500KHZ.svd ';  %189 13
fileName = 'Scan_time_Allu_500KHZ.svd ';

filename = sprintf('%s%s',pathName,fileName);
eval(sprintf('cd %s',pathName));
%eval(sprintf('cd %s',pathcurrent));

%--------------------------------------------------------------------------------------------------------------------
%                       Get data points from the svd or pvd file
%--------------------------------------------------------------------------------------------------------------

% we use the fuction GetPointData, for the détails see Anex s.1
domainname='Time';
channelname='Vib X';
% signalname='Displacement';
signalname='Velocity';
displayname='Samples';
point=0;
frame=0;
% fftAquis = GetFFTAcqProps(filename);
[t,vx,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point, frame);
channelname='Vib Y';
[t,vy,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point, frame);
channelname='Vib Z';
[t,vz,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point, frame);

%Remark: for a 3D acquisition on nx*ny points on nt time points, vx is organized as follow:
%   - vx is a 2D matrix: nx*ny*nz lines and nt columns
%   - nx,ny are required for reshaping data
figure(1)
subplot(2,1,1)
hold on
m=20;
plot(t,vx(m,:),t,vy(m,:),'r',t,vz(m,:),'k');
title('vx-vy-vz')

subplot(2,1,2)
imagesc(vy)

%---------------------------------------------------------------------------------------------------
%                           Matrix of Index XYZCoordinates 
%-----------------------------------------------------------------------------------------
%
%We use the function GetXYZCoordinates(filename, point) see Anex s.2
IndexCordinat = GetXYZCoordinates(filename, point);
%%
figure(2)
subplot(2,1,1)
%new from 24/6
XYZ= GetXYZCoordinates(filename, point);
 surf(reshape(XYZ(:,1),189,13),reshape(XYZ(:,2),189,13),reshape(XYZ(:,3),189,13),reshape(vx(:,1900),189,13))
subplot(2,1,2)
for i=1:200, plot3(reshape(XYZ(:,1),189,13),reshape(XYZ(:,2),189,13),reshape(vx(:,i),189,13),'.'),i,pause,end
nx=189
ny=13
nt=5000
for i=1:200, imagesc(reshape(vx(:,i),nx,ny)),i,pause,end
save Rayleigh3D
else
    load Rayleigh3D.mat
end
%%
%Caution: z on propa axis, x normal to surf.... for ALU

vvx=flipud(vz);
vvy=flipud(vy);
vvz=flipud(vx);
[nx nt]=size(vx);
%goto complex?
% cv=squeeze(vvx(:,5,1:600))+sqrt(-1)*squeeze(vvz(:,5,1:600));
%Caution: z on propa axis, x normal to surf.... for ALU
cv=sqrt(-1)*squeeze(vvx(:,1:600))+squeeze(vvz(:,1:600));
%cvh=hilbert(squeeze(vvx(:,1:600)))+sqrt(-1)*hilbert(squeeze(vvz(:,1:600)));

% vvx=reshape(vx,nx,ny,nt);
% vvy=reshape(vy,nx,ny,nt);
% vvz=reshape(vz,nx,ny,nt);

dx=abs(min(diff(XYZ(:,3))));
x=(0:nx-1)*dx;
% dy=abs(min(diff(XYZ(:,2))));
% y=(0:ny-1)*dy;
% dz=abs(min(diff(XYZ(:,3))));
dt=min(diff(t))
f=((0:nt-1)/nt-.5)/dt;

%Compare 3 directions on all FOV:
figure(3)
mmax=max([vvx(:);vvy(:);vvz(:)])/10;
subplot(231)
    imagesc(t,x,vvx,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('x vibration')
subplot(232)
    imagesc(t,x,vvy,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('y vibration')
subplot(233)
    imagesc(t,x,vvz,[0 mmax])
    xlabel('Time (s)')
    ylabel('x direction')
    title('z vibration')
    colorbar
  load config_file
  load signal_600_1.mat
    tmp=resample(T,1,10);
    nt=length(tmp)-1;
    dt=median(diff(tmp));
    t0=T(1);
    t=t0+(0:nt-1)*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx=600;
    ny=1;
    dx=0.2e-3;
    dt=median(diff(t));
    x=(0:nx-1)*dx;
    k=(((0:nx-1)/nx)-.5)/dx;
    f=(((0:nt-1)/nt)-.5)/dt;
subplot(234)
    ix=imagesc(f,-k,fftshift(abs(fft2(vvx)))),hold,plot(f,f./2930),hold;
    xlabel(' f ')
    ylabel('k ')
    title('Spectre 2D et comparaison')
    colorbar
subplot(235)
    iy=imagesc(f,-k,fftshift(abs(fft2(vvy)))),hold,plot(f,f./2930),hold;
    xlabel(' f ')
    ylabel('k ')
    title('Spectre 2D et comparaison')
    colorbar
subplot(236)
    iz=imagesc(f,-k,fftshift(abs(fft2(vvz)))),hold,plot(f,f./2930),hold;
    xlabel(' f ')
    ylabel('k ')
    title('Spectre 2D et comparaison')
    colorbar
        
figure(4)
sigma = 1.6;
subplot(2,2,1)
RR =fftshift(abs(fft2(vvx)));% fft
imagesc(RR);
% se = strel('disk',3);
% RR=imclose(RR,se);
RR = wiener2(RR,[5 10]); %filtre wiener2D
subplot(2,2,2)
imagesc(RR);
thresh=[0.01,0.17];
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
subplot(2,2,3)
f=imdilate(f,ones(3));
imagesc(f)
title('Bord ')
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',25:0.01:80); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
%%
figure(5)
lines = houghlines(f,T,R,P,'FillGap',2,'Minlength',20); % lignes
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    len = norm(lines(k).point1 - lines(k).point2);
    Len(k)=len;
    if(len>45 & len<200)
        plot(xy(:,1),-xy(:,2),'LineWidth',1,'color','r');
        
        axis([0 2500,-437 0])
        
        set(gca,'xtick',[],'ytick',[]);
        %axis off;
        hold on;
    end
end



imgx = fftshift(fft2(vvx));
saveas(gcf,'filt','png');
fil = imread('filt.png');
ci=imresize(fil,[437,2500]);
ci = rgb2gray(ci);
level=graythresh(ci)
ci = im2bw(ci,level);
imgx(ci==1)=0;
figure(6)
subplot(121)
imshow(imgx);
%%
bwn = imgx(1:ceil(size(imgx,1)/2),size(imgx,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
subplot(122)

%%
for i=1200:1300
    fill(:,i)=0;
end
imagesc(abs(fill));
imageo = abs(ifft2(ifftshift(fill)));
figure(7)
imagesc(imageo);
% for i=1:437
%     for j=1:2500
%         if c(i,j)==0
%             n=n+1;
%             i+1;
%         end
%     end
% end
 %%       
%       figure(5)
%
% figure(4)
% subplot(221)
%     plot(t,vvx(200,:))
%     xlabel('Time (s)')
%     ylabel('Vibration Vx')
%     title('Vx(t) for x=44.6mm')
%     axis([0 5e-5 -5e-3 5e-3])
% subplot(222)
%     plot(t,vvy(200,:))
%     xlabel('Time (s)')
%     ylabel('Vibration Vy')
%     title('Vy(t) for x=44.6mm')
%     axis([0 5e-5 -5e-3 5e-3])
% subplot(223)
%     plot(t,vvz(200,:))
%     xlabel('Time (s)')
%     ylabel('Vibration Vz')
%     title('Vz(t) for x=44.6mm')
%     axis([0 5e-5 -5e-3 5e-3])
% 
% %All
% subplot(224)
%     plot(t,vvx(200,:),t,vvy(200,:),t,vvz(200,:))
%     xlabel('Time (s)')
%     ylabel('Vibrations Vx, Vy, Vz')
%     title('V for x=44.6mm')
%     axis([0 5e-5 -5e-3 5e-3])
% 
% %Ellipticity:
% %theoretical ellipso H/V    
% qt=(2765/3388)^2;
% hvt=2*(sqrt(1-qt)./(2-qt)) %.866
% 
% %%
% figure(5)
% 
% int=120:190;
% subplot(211)
%     plot(t(int),vvx(200,int),t(int),vvz(200,int))
%     xlabel('Time (s)')
%     ylabel('Vibrations Vx, Vz')
%     title('V for x=44.6mm')
%     %axis([0 5e-5 -5e-3 5e-3])
% subplot(212)
% mmax=max([vvx(200,int),vvz(200,int)]');
% mvvx=vvx(200,int)-mean(vvx(200,int));
% mvvz=vvz(200,int)-mean(vvz(200,int));
% mmax=max([mvvx,mvvz]');
% 
%     [zet]=elliptic_plot(mmax*hvt,mmax,0,0,0);
%     [ze]=elliptic_plot(mmax*.5,mmax,0,0,0);
%     plot(vvx(200,int),vvz(200,int),real(zet),imag(zet),real(ze),imag(ze))
%     xlabel('Vx')
%     ylabel('Vz')
%     title('Ellipsometry for x=44.6mm')
%     axis([-5e-3 5e-3 -5e-3 5e-3])
%     axis equal
% %theoretical ellipso H/V    
% q=(2765/3388)^2;
% hv=2*(sqrt(1-q)./(2-q))
% 
