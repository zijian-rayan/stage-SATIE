
% Checking Polytec import Matlab functions for svd files:
%(1)Prior use, install Polytec libraries:
%   - Polytec_Common_Runtime_171_64Bit_(1)
%   - Polytec_FileAccess_471_(2)


clear all;
close all;
clc;
importfiles=0;
%if importfiles==1,
    
% % opengl('save','software');
% global pathcurrent;
% pathcurrent = pwd;
% %pathName = 'C:\Users\BILY_HOLLOW\Desktop\Stage_2015\DATA\PSV_3D\Marbre\Format_svd\';
% pathName = 'D:\ÊµÑéÊÒÊµÏ°\Signaux1D_Rayleigh\';
% 
% 
% %fileName = 'Scan_time_Marbre_g1_Alt_500KHZ.svd';
% %fileName = 'Scan_time_Marbre_g2_Alt_500KHZ.svd';
% %fileName = 'Scan_time_MarbreI1_500KHZ.svd ';  %189 13
% fileName = 'Scan_time_Allu_500KHZ.svd ';
% 
% filename = sprintf('%s%s',pathName,fileName);
% eval(sprintf('cd %s',pathName));
% %eval(sprintf('cd %s',pathcurrent));
% 
% %--------------------------------------------------------------------------------------------------------------------
% %                       Get data points from the svd or pvd file
% %--------------------------------------------------------------------------------------------------------------
% 
% % we use the fuction GetPointData, for the détails see Anex s.1
% domainname='Time';
% channelname='Vib X';
% % signalname='Displacement';
% signalname='Velocity';
% displayname='Samples';
% point=0;
% frame=0;
% % fftAquis = GetFFTAcqProps(filename);
% [t,vx,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point, frame);
% channelname='Vib Y';
% [t,vy,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point, frame);
% channelname='Vib Z';
% [t,vz,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point, frame);
% 
% %Remark: for a 3D acquisition on nx*ny points on nt time points, vx is organized as follow:
% %   - vx is a 2D matrix: nx*ny*nz lines and nt columns
% %   - nx,ny are required for reshaping data
% figure(1)
% subplot(2,1,1)
% hold on
% m=20;
% plot(t,vx(m,:),t,vy(m,:),'r',t,vz(m,:),'k');
% title('vx-vy-vz')
% 
% subplot(2,1,2)
% imagesc(vy)
% 
% %---------------------------------------------------------------------------------------------------
% %                           Matrix of Index XYZCoordinates 
% %-----------------------------------------------------------------------------------------
% %
% %We use the function GetXYZCoordinates(filename, point) see Anex s.2
% IndexCordinat = GetXYZCoordinates(filename, point);
% %%
% figure(2)
% subplot(2,1,1)
% %new from 24/6
% XYZ= GetXYZCoordinates(filename, point);
%  surf(reshape(XYZ(:,1),189,13),reshape(XYZ(:,2),189,13),reshape(XYZ(:,3),189,13),reshape(vx(:,1900),189,13))
% subplot(2,1,2)
% for i=1:position, plot3(reshape(XYZ(:,1),189,13),reshape(XYZ(:,2),189,13),reshape(vx(:,i),189,13),'.'),i,pause,end
% nx=189
% ny=13
% nt=5000
% for i=1:position, imagesc(reshape(vx(:,i),nx,ny)),i,pause,end
% save Rayleigh3D
% else
%     load Rayleigh3D.mat
% end
 load Rayleigh3D.mat
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
% figure(3)
% mmax=max([vvx(:);vvy(:);vvz(:)])/10;
% subplot(231)
%     imagesc(t,x,vvx,[0 mmax])
%     colorbar
%     xlabel('Time (s)')
%     ylabel('x direction')
%     title('x vibration')
% subplot(232)
%     imagesc(t,x,vvy,[0 mmax])
%     colorbar
%     xlabel('Time (s)')
%     ylabel('x direction')
%     title('y vibration')
% subplot(233)
%     imagesc(t,x,vvz,[0 mmax])
%     xlabel('Time (s)')
%     ylabel('x direction')
%     title('z vibration')
%     colorbar
%   load config_file
%   load signal_600_1.mat
%     tmp=resample(T,1,10);
%     nt=length(tmp)-1;
%     dt=median(diff(tmp));
%     t0=T(1);
%     t=t0+(0:nt-1)*dt;
% %%%%%%%%%%%%%%%%%%%%%%%%%%
%     nx=600;
%     ny=1;
%     dx=0.2e-3;
%     dt=median(diff(t));
%     x=(0:nx-1)*dx;
%     k=(((0:nx-1)/nx)-.5)/dx;
%     f=(((0:nt-1)/nt)-.5)/dt;
% subplot(234)
%     ix=imagesc(f,-k,fftshift(abs(fft2(vvx)))),hold,plot(f,f./2930),hold;
%     xlabel(' f ')
%     ylabel('k ')
%     title('Spectre 2D et comparaison')
%     colorbar
% subplot(235)
%     iy=imagesc(f,-k,fftshift(abs(fft2(vvy)))),hold,plot(f,f./2930),hold;
%     xlabel(' f ')
%     ylabel('k ')
%     title('Spectre 2D et comparaison')
%     colorbar
% subplot(236)
%     iz=imagesc(f,-k,fftshift(abs(fft2(vvz)))),hold,plot(f,f./2930),hold;
%     xlabel(' f ')
%     ylabel('k ')
%     title('Spectre 2D et comparaison')
%     colorbar
        
figure(4)
sigma = 1.6;
subplot(131)
RR =fftshift(abs(fft2(vvx)));% fft
imagesc(RR);
% se = strel('disk',3);
% RR=imclose(RR,se);
RR = wiener2(RR,[5 10]); %filtre wiener2D
subplot(132)
imagesc(RR);
thresh=[0.11,0.37];
% thresh=[0.01,0.17];
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
subplot(133)
%f=imdilate(f,ones(3));
imagesc(f)
title('Bord ')
%[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);

%%
% figure(5)
% lines = houghlines(f,T,R,P,'FillGap',5,'Minlength',20); % lignes
% for k = 1:length(lines)
%     xy = [lines(k).point1; lines(k).point2];
%     len = norm(lines(k).point1 - lines(k).point2);
%     Len(k)=len;
%     if(len>20 & len<1000)
%         plot(xy(:,1),-xy(:,2),'LineWidth',1,'color','b');
%         
%         axis([0 2500,-437 0])
%         
%         box off
%         ax2 = axes('Position',get(gca,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
%             set(ax2,'YTick', []);
%             set(ax2,'XTick', []);
%             box on
%         %set(gca,'xtick',[],'ytick',[]);
%         %axis off;
%         hold on;
%     end
% end
% saveas(gcf,'filt','png');
%%
% figure(6)
% fil = imread('filt.png');
% ci=imresize(fil,[437,2500]);
% ci = rgb2gray(ci);
% level=graythresh(ci)
% ci = im2bw(ci,level);
% se = strel('disk',1);
% ci=imdilate(abs(ci),se);
% bwn = ci(1:ceil(size(ci,1)/2),size(ci,2)/2+1:end);
% bwl = fliplr(bwn);
% bwd = flipud(bwn);
% bwld = fliplr(bwd);
% fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
% imaxx = fftshift((fft2(vvx)));
% imaxx(fill==1)=0;
% subplot(321)
% imagesc(abs(imaxx));
% title('Image  FFT');
% colorbar;
% subplot(322)
% %6-28%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imageo = abs(ifft2(ifftshift(imaxx)));
% imagesc(imageo);
% title('Image filtr¨¦e');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fil = imread('filt.png');
% ci=imresize(fil,[437,2500]);
% ci = rgb2gray(ci);
% level=graythresh(ci)
% ci = im2bw(ci,level);
% se = strel('disk',3);
% ci=imdilate(abs(ci),se);
% bwn = ci(1:ceil(size(ci,1)/2),size(ci,2)/2+1:end);
% bwl = fliplr(bwn);
% bwd = flipud(bwn);
% bwld = fliplr(bwd);
% fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
% imaxx = fftshift((fft2(vvx)));
% imaxx(fill==1)=0;
% subplot(323)
% imagesc(abs(imaxx));
% title('Image  FFT');
% colorbar;
% subplot(324)
% imageo = abs(ifft2(ifftshift(imaxx)));
% imagesc(imageo);
% title('Image filtr¨¦e');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fil = imread('filt.png');
% ci=imresize(fil,[437,2500]);
% ci = rgb2gray(ci);
% level=graythresh(ci)
% ci = im2bw(ci,level);
 se = strel('disk',5);
% ci=imdilate(abs(ci),se);
% bwn = ci(1:ceil(size(ci,1)/2),size(ci,2)/2+1:end);
% bwl = fliplr(bwn);
% bwd = flipud(bwn);
% bwld = fliplr(bwd);
% fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
% imaxx = fftshift((fft2(vvx)));
% imaxx(fill==1)=0;
% subplot(325)
% imagesc(abs(imaxx));
% title('Image  FFT');
% colorbar;
% subplot(326)
% imageo = abs(ifft2(ifftshift(imaxx)));
% imagesc(imageo);
% title('Image filtr¨¦e');
%%
figure(8)
subplot(221)
imagesc(f);
colorbar
title('Fronti¨¨re')
f=imdilate(abs(f),se);
f=imerode(f,se);
    f(:,1225:1275)=0;
    f(1:50,:)=0;
    f(387:437,:)=0;

subplot(222)
imagesc(f);
title('Filtre trait¨¦')
colorbar
imaxx = fftshift((fft2(vvx)));
imaxx(f==0)=0;
subplot(223)
imagesc(abs(imaxx));
title('Image de fr¨¦quence filtr¨¦e')
subplot(224)
imageo = real(ifft2(ifftshift(imaxx)));
imagesc(imageo);
title('Image en domaine temporel')
%%
figure(9)
bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
subplot(221)
imagesc(fill);
colorbar
title('image fill')
subplot(222)
sigma = 1;
gausFilter = fspecial('gaussian', [5,5], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imagesc(gaus);
colorbar
title('gaussian')
subplot(223)
imaxx = fftshift((fft2(vvy)));
imaxx=imaxx.*gaus;
imagesc(abs(imaxx))
colorbar
title('imaxx gaussian')
subplot(224)
imageo2 = real(ifft2(ifftshift(imaxx)));
imagesc(imageo2);
colorbar
title('Image en domaine temporel')
figure(10)
imaxor = imageo2-imageo;
imagesc(imaxor);
colorbar
title('Image en domaine temporel')
 s=median(imageo2);
mean(s)

 %%       
%    
% 
position = 300;
figure(9)
subplot(221)
    plot(t,vvx(position,:))
    xlabel('Time (s)')
    ylabel('Vibration Vx')
    title('Vx(t) for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])
subplot(222)
    plot(t,vvy(position,:))
    xlabel('Time (s)')
    ylabel('Vibration Vy')
    title('Vy(t) for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])
subplot(223)
    plot(t,vvz(position,:))
    xlabel('Time (s)')
    ylabel('Vibration Vz')
    title('Vz(t) for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])

%All
subplot(224)
    plot(t,vvx(position,:),t,vvy(position,:),t,vvz(position,:))
    xlabel('Time (s)')
    ylabel('Vibrations Vx, Vy, Vz')
    title('V for x=44.6mm')
    axis([0 5e-5 -5e-3 5e-3])

%Ellipticity:
%theoretical ellipso H/V    
qt=(2765/3388)^2;
hvt=2*(sqrt(1-qt)./(2-qt)) %.866


figure(10)

int=120:190;
subplot(211)
    plot(t(int),vvx(position,int),t(int),vvz(position,int))
    xlabel('Time (s)')
    ylabel('Vibrations Vx, Vz')
    title('V for x=44.6mm')
    %axis([0 5e-5 -5e-3 5e-3])
subplot(212)
mmax=max([vvx(position,int),vvz(position,int)]');
mvvx=vvx(position,int)-mean(vvx(position,int));
mvvz=vvz(position,int)-mean(vvz(position,int));
mmax=max([mvvx,mvvz]');

    [zet]=elliptic_plot(mmax*hvt,mmax,0,0,0);
    [ze]=elliptic_plot(mmax*.5,mmax,0,0,0);
    plot(vvx(position,int),vvz(position,int))
    xlabel('Vx')
    ylabel('Vz')
    title('Ellipsometry for x=44.6mm')
    axis([-5e-3 5e-3 -5e-3 5e-3])
    axis equal
%theoretical ellipso H/V    
q=(2765/3388)^2;
hv=2*(sqrt(1-q)./(2-q))
% 
