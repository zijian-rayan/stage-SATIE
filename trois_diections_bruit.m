%direction y bruit
% Ce fichier contient les effets de filtrage des signaux dans trois directions dans les domaines temps et fr¨¦quence.
% Et s'ins¨¦rant dans une ellipse
clear all;
close all;
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
%  se = strel('disk',5);
%  filt_y1=imdilate(abs(filt_y1),se);
% filt_y1=imerode(filt_y1,se);


% for i=1:437
%         B = sort(unique(vvy(i,:)), 'descend');
%         for j=1:2500
%             if abs(filt_y2(i,j))<B(30)
%                 filt_y2(i,j)=0;
%             else
%                 if filt_y2(i,j)>0
%                     filt_y2(i,j)=filt_y2(i,j)-B(30);
%                 else
%                     filt_y2(i,j)=filt_y2(i,j)+B(30);
%                 end
%             end
%         end
% end
% for i=1:437
%             B = sort(unique(vvy(i,:)), 'descend');
%         for j=1:2500
%             if abs(filt_y3(i,j))<B(20)
%                 filt_y3(i,j)=0;
%             else
%                 if filt_y3(i,j)>0
%                     filt_y3(i,j)=filt_y3(i,j)-B(20);
%                 else
%                     filt_y3(i,j)=filt_y3(i,j)+B(20);
%                 end
%             end
%         end
%         for j=1:2500
%             if abs(filt_y3(i,j))<ave(i)
%                 filt_y3(i,j)=0;
%             else
%                 if filt_y3(i,j)>0
%                     filt_y3(i,j)=filt_y3(i,j)-ave(i);
%                 else
%                     filt_y3(i,j)=filt_y3(i,j)+ave(i);
%                 end
%             end
%         end
%end
%
%     figure(2)
%     nx=600;
%     ny=1;
%     dx=0.2e-3;
%     dt=median(diff(t));
%     x=(0:nx-1)*dx;
%     k=(((0:nx-1)/nx)-.5)/dx;
%     f=(((0:nt-1)/nt)-.5)/dt;
% subplot(221)
%     iy=imagesc(f,-k,fftshift(abs(fft2(vvy))));
%     xlabel(' f ')
%     ylabel('k ')
%     title('Spectre 2D et comparaison')
%     colorbar
%     figure(2)
% subplot(222)
%     iy=imagesc(f,-k,fftshift(abs(fft2(filt_y1))));
%     xlabel(' f ')
%     ylabel('k ')
%     title('Spectre 2D et comparaison')
%     colorbar
%         figure(2)
% subplot(223)
%     iy=imagesc(f,-k,fftshift(abs(fft2(filt_y2))));
%     xlabel(' f ')
%     ylabel('k ')
%     title('Spectre 2D et comparaison')
%     colorbar
%         figure(2)
% subplot(224)
%     iy=imagesc(f,-k,fftshift(abs(fft2(filt_y3))));
%     xlabel(' f ')
%     ylabel('k ')
%     title('Spectre 2D et comparaison')
%     colorbar
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
%subplot(131)
%imagesc(RR)
%title('y direction')
%subplot(132)
%imagesc(abs(imaxx))
bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
%fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
fill = f;
sigma = 1;
    fill(:,1180:1320)=0;
    fill(1:50,:)=0;
    fill(387:437,:)=0;
gausFilter = fspecial('gaussian', [3,3], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(filt_y1)));
imaxx=imaxx.*gaus;
imageoy = real(ifft2(ifftshift(imaxx)));
%subplot(133)
%imagesc(abs(imaxx))
%%
figure(2)
subplot(345)
    imagesc(t,x,vvy,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Image originale en direction Y')
subplot(346)
    imagesc(t,x,filt_y1,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Filtrage dans le domaine temporel dans la direction Y')
subplot(347)
    imagesc(t,x,imageoy,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Filtrage du domaine fr¨¦quentiel dans la direction Y')
subplot(348)
    mesh(imageoy)
   
s=var(imageoy);
mean_var_y = mean(s)
mean_mean_y = mean(mean(imageoy))
mean_max_y = mean(max(imageoy))
mean_median_y = mean(median(imageoy))
%%
%figure(3)
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
% subplot(131)
% imagesc(RR)
% title('z direction')
% subplot(132)
% imagesc(abs(imaxx))
bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
%fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
fill = f;
sigma = 1;
    fill(:,1180:1320)=0;
    fill(1:50,:)=0;
    fill(387:437,:)=0;
gausFilter = fspecial('gaussian', [3,3], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(filt_z)));
imaxx=imaxx.*gaus;
imageoz = real(ifft2(ifftshift(imaxx)));
%subplot(133)
%imagesc(abs(imaxx))
%%
figure(2)
subplot(3,4,9)
    imagesc(t,x,vvz,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Image originale en direction Z')
subplot(3,4,10)
    imagesc(t,x,filt_z,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Filtrage dans le domaine temporel dans la direction Z')
subplot(3,4,11)
    imagesc(t,x,imageoz,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Filtrage du domaine fr¨¦quentiel dans la direction Z')
subplot(3,4,12)
    mesh(imageoz)
s=var(imageoz);
mean_var_z = mean(s)
mean_mean_z = mean(mean(imageoz))
mean_max_z = mean(max(imageoz))
mean_median_z = mean(median(imageoz))
%%
%figure(5)
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
% subplot(131)
% imagesc(RR)
% title('x direction')
% subplot(132)
% imagesc(abs(imaxx))
bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
%fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
fill = f;
sigma = 1;
    fill(:,1180:1320)=0;
    fill(1:50,:)=0;
    fill(387:437,:)=0;
gausFilter = fspecial('gaussian', [3,3], sigma);
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(filt_x)));
imaxx=imaxx.*gaus;
imageox = real(ifft2(ifftshift(imaxx)));
%subplot(133)
%imagesc(abs(imaxx))
%%
figure(2)
subplot(341)
    imagesc(t,x,vvx,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Image originale en direction X')
subplot(3 ,4, 2)
    imagesc(t,x,filt_x,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Filtrage du domaine fr¨¦quentiel dans la direction X')
subplot(3 ,4 ,3)
    imagesc(t,x,imageox,[0 mmax])
    colorbar
    xlabel('Time (s)')
    ylabel('x direction')
    title('Filtrage du domaine fr¨¦quentiel dans la direction X')
subplot(3, 4, 4)
    mesh(imageox)
s=var(imageox);
mean_var_x = mean(s)
mean_mean_x = mean(mean(imageox))
mean_max_x = mean(max(imageox))
mean_median_x = mean(median(imageox))
% figure
% imagesc(abs(imaxx));

%%
% %%%%%%%%%%%%%%%%%%%%% y2
% sigma = 1.6;
% RR =fftshift(abs(fft2(filt_y1)));% fft
% RR = wiener2(RR,[5 10]); %filtre wiener2D
% thresh=[0.11,0.121];
% % thresh=[0.01,0.17];
% f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
% %[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
% [H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
% P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
% se = strel('disk',3);
% se2 = strel('disk',5);
% f=imdilate(abs(f),se);
% f=imerode(f,se2);
%     f(:,1180:1320)=0;
%     f(1:50,:)=0;
%     f(387:437,:)=0;
% imaxx = fftshift((fft2(filt_y2)));
% imaxx(f==0)=0;
% imageo = real(ifft2(ifftshift(imaxx)));
% 
% bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
% bwl = fliplr(bwn);
% bwd = flipud(bwn);
% bwld = fliplr(bwd);
% fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
% sigma = 1;
%     fill(:,1180:1320)=0;
%     fill(1:50,:)=0;
%     fill(387:437,:)=0;
% gausFilter = fspecial('gaussian', [5,5], sigma);
% gaus= imfilter(fill, gausFilter, 'replicate');
% imaxx = fftshift((fft2(filt_y2)));
% imaxx=imaxx.*gaus;
% imageo2 = real(ifft2(ifftshift(imaxx)));
% imaxor = imageo2-imageo;
% subplot(334)
%     imagesc(t,x,filt_y2,[0 mmax])
%     colorbar
%     xlabel('Time (s)')
%     ylabel('x direction')
%     title('y vibration')
% subplot(335)
%     imagesc(t,x,imageo2,[0 mmax])
%     colorbar
%     xlabel('Time (s)')
%     ylabel('x direction')
%     title('y vibration')
% subplot(336)
%     mesh(imageo2)
% s=var(imageo2);
% mean_y = mean(s)
% 
% %%
% %%%%%%%%%%%%%%%%%y3
% 
% sigma = 1.6;
% RR =fftshift(abs(fft2(filt_y1)));% fft
% RR = wiener2(RR,[5 10]); %filtre wiener2D
% thresh=[0.11,0.121];
% % thresh=[0.01,0.17];
% f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
% %[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10); %tr. hough
% [H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
% P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
% se = strel('disk',3);
% se2 = strel('disk',5);
% f=imdilate(abs(f),se);
% f=imerode(f,se2);
%     f(:,1180:1320)=0;
%     f(1:50,:)=0;
%     f(387:437,:)=0;
% imaxx = fftshift((fft2(filt_y3)));
% imaxx(f==0)=0;
% imageo = real(ifft2(ifftshift(imaxx)));
% 
% bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);
% bwl = fliplr(bwn);
% bwd = flipud(bwn);
% bwld = fliplr(bwd);
% fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];
% sigma = 1;
%     fill(:,1180:1320)=0;
%     fill(1:50,:)=0;
%     fill(387:437,:)=0;
% gausFilter = fspecial('gaussian', [5,5], sigma);
% gaus= imfilter(fill, gausFilter, 'replicate');
% imaxx = fftshift((fft2(filt_y3)));
% imaxx=imaxx.*gaus;
% imageo2 = real(ifft2(ifftshift(imaxx)));
% imaxor = imageo2-imageo;
% subplot(337)
%     imagesc(t,x,filt_y3,[0 mmax])
%     colorbar
%     xlabel('Time (s)')
%     ylabel('x direction')
%     title('y vibration')
% subplot(338)
%     imagesc(t,x,imageo2,[0 mmax])
%     colorbar
%     xlabel('Time (s)')
%     ylabel('x direction')
%     title('y vibration')
% subplot(339)
%     mesh(imageo2)
% s=var(imageo2);
% mean_y = mean(s)
% 

%%
% figure
% clf;
% for i = 1:437
%     subplot(221)
%     plot(imageo2(:,i),'b');
%     title(i)
%     text(30,30,'Scannez chaque ligne')
%     %hold on;
%     axis([0 437 -0.0006 0.0006]);
%     subplot(222)
%     plot(imageo2(i,:),'r');
%     %hold on;
%     axis([0 2500 -0.0005 0.0005]);
%     pause(0.001);
%     title(i)
%     text(300,30,'Scannez chaque colonne')
%         subplot(223)
%     plot(imageo2(:,i),'b');
%     title(i)
%     text(30,30,'Scannez chaque ligne')
%     hold on;
%     axis([0 437 -0.0006 0.0006]);
%     subplot(224)
%     plot(imageo2(i,:),'r');
%     hold on;
%     axis([0 2500 -0.0005 0.0005]);
%     pause(0.001);
%     title(i)
%     text(300,30,'Scannez chaque colonne')
% end



%%
figure(4)
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

%%
figure(6)

int=120:190;
subplot(131)
mmax=max([imageox(200,int),imageoz(200,int)]');
mvvx=imageox(200,int)-mean(imageox(200,int));
mvvz=imageoz(200,int)-mean(imageoz(200,int));
mmax=max([mvvx,mvvz]');

    [zet]=elliptic_plot(mmax*hvt,mmax,0,0,0);
    [ze]=elliptic_plot(mmax*.5,mmax,0,0,0);
    plot(imageox(200,int),imageoz(200,int))
    xlabel('Vx')
    ylabel('Vz')
    title('Ellipsometry X-Z for x=44.6mm')
    axis([-1e-3 1e-3 -1e-3 1e-3])
    axis equal
%theoretical ellipso H/V    
q=(2765/3388)^2;
hv=2*(sqrt(1-q)./(2-q))

subplot(132)
mmax=max([imageox(200,int),imageoy(200,int)]');
mvvx=imageox(200,int)-mean(imageox(200,int));
mvvz=imageoy(200,int)-mean(imageoy(200,int));
mmax=max([mvvx,mvvz]');

    [zet]=elliptic_plot(mmax*hvt,mmax,0,0,0);
    [ze]=elliptic_plot(mmax*.5,mmax,0,0,0);
    plot(imageox(200,int),imageoy(200,int))
    xlabel('Vx')
    ylabel('Vz')
    title('Ellipsometry X-Y for x=44.6mm')
    axis([-1e-3 1e-3 -1e-3 1e-3])
    axis equal

subplot(133)
mmax=max([imageoz(200,int),imageoy(200,int)]');
mvvx=imageoz(200,int)-mean(imageoz(200,int));
mvvz=imageoy(200,int)-mean(imageoy(200,int));
mmax=max([mvvx,mvvz]');

    [zet]=elliptic_plot(mmax*hvt,mmax,0,0,0);
    [ze]=elliptic_plot(mmax*.5,mmax,0,0,0);
    plot(imageoz(200,int),imageoy(200,int))
    xlabel('Vx')
    ylabel('Vz')
    title('Ellipsometry Y-Z for x=44.6mm')
    axis([-1e-3 1e-3 -1e-3 1e-3])
    axis equal



