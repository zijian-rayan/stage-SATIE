%% prony

clear all;
close all;
clc;
load Rayleigh3D.mat
[nx nt]=size(vx);
dx=abs(min(diff(XYZ(:,3))));%D¨¦finir l'intervalle d'¨¦chantillonnage
x=(0:nx-1)*dx;
vvy=flipud(vy);%Lire le signal
vvx=flipud(vz);
vvz=flipud(vx);
mmax=max([vvy(:)])/10;
filt_z=vvz;
filt_y1=vvy;
filt_x=vvx;
for i=1:437%Filtrage dans le domaine temporel
    B = sort(unique(vvx(i,:)), 'descend');
    M = mean(vvx(i,:));
    for j=1:2500
        filt_x(i,j)=filt_x(i,j)-M;
        if abs(filt_x(i,j))<B(20)
            filt_x(i,j)=0;
        else
            if filt_x(i,j)>0
                filt_x(i,j)=filt_x(i,j)-B(20);  %Compensation du signal
            else
                filt_x(i,j)=filt_x(i,j)+B(20);
            end
        end
    end
end
sigma = 1.6;
RR =fftshift(abs(fft2(filt_x)));% fft
RR = wiener2(RR,[5 10]); %filtre wiener2D
thresh=[0.11,0.37];%Seuil, cette valeur est diff¨¦rente dans trois directions
f=edge(double(RR),'canny',thresh,sigma); % obtenir le bord
[H,T,R]=hough(f,'Theta',5:0.01:45); %tr. hough
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);
se2 = strel('disk',5);%Cr¨¦er un masque
f=imdilate(abs(f),se2);%Remplissage de trous
f=imerode(f,se2);
f(:,1180:1320)=0;%?liminer le bruit de basse fr¨¦quence
f(1:50,:)=0;
f(387:437,:)=0;
imaxx = fftshift((fft2(filt_x)));
imaxx(f==0)=0;%Cr¨¦er le filtre de domaine de fr¨¦quence
imageo = real(ifft2(ifftshift(imaxx)));
bwn = f(1:ceil(size(f,1)/2),size(f,2)/2+1:end);%Retournez le filtre qui g¨¦n¨¨re quatre quadrants
bwl = fliplr(bwn);
bwd = flipud(bwn);
bwld = fliplr(bwd);
fill = [bwl(1:end-1,:),bwn(1:end-1,:);bwld,bwd];%Filtre complet
sigma = 1;
fill(:,1180:1320)=0;%?liminer le bruit de basse fr¨¦quence
fill(1:50,:)=0;
fill(387:437,:)=0;
gausFilter = fspecial('gaussian', [3,3], sigma);%Lissage des bords
gaus= imfilter(fill, gausFilter, 'replicate');
imaxx = fftshift((fft2(filt_x)));
imaxx=imaxx.*gaus;%Filtrer
imageox = real(ifft2(ifftshift(imaxx)));%R¨¦sultat du filtrage dans le domaine temporel

%%
%$$$$$ prony  $$$$


Fs=10e6
dt=1/Fs
dx=abs(min(diff(XYZ(:,3))));%D¨¦finir l'intervalle d'¨¦chantillonnage
Fx=1/dx
s = real(imageox);%Chargement du signal
SS=fft(s,[],2);
freq=((0:nt-1)/(nt)-.5)*2*pi;
fff=((0:nt-1)/(nt)-.5)*Fs;
kkk=fliplr((0:nx-1)/(nx)-.5)*Fx;
clear rr pp kk aa bb%M¨¦moire vide
for i=1:nt%M¨¦thode prony
	[b,a]=prony(SS(:,i),3,4);
	aa(i,1:length(a))=a;%Mol¨¦cule et d¨¦nominateur
	bb(i,1:length(b))=b;
   [r,p,k]=residue(b,a);%D¨¦composer en p?les et r¨¦sidus
	rr(i,1:length(r))=r';
	pp(i,1:length(p))=p';
   kk(i,1:length(k))=k';
end
VPest=log(abs(pp(:,1)))*Fx;%Chargement de poteau
intex = find(VPest>-100);
fff=fftshift(fff);
figure(1)
plot((fff(1,intex)),VPest(intex,1),'b.');hold on;%Pole print
xlabel('Fr¨¦quence (Hz)')
ylabel('Pole')
title('Pole')
figure(2)
fff=fftshift(fff);
VPest=abs(2*pi*fff'./(fftshift(angle(pp(:,2))*Fx)+1e-6));
semilogy(fff(1,1050:1175),VPest(1050:1175,1),'b.',fff(1,1325:1450),VPest(1325:1450,1),'b.');hold on;%Vitesse de phase
xlabel('Fr¨¦quence (Hz)')
ylabel('Vitesse de phase (m/s)')
title('imageox')
%%
s = real(vvx);
SS=fft(s,[],2);
freq=((0:nt-1)/(nt)-.5)*2*pi;
fff=((0:nt-1)/(nt)-.5)*Fs;
kkk=fliplr((0:nx-1)/(nx)-.5)*Fx;

clear rr pp kk aa bb
for i=1:nt
	[b,a]=prony(SS(:,i),3,4);
	aa(i,1:length(a))=a;
	bb(i,1:length(b))=b;
   [r,p,k]=residue(b,a);
	rr(i,1:length(r))=r';
	pp(i,1:length(p))=p';
   kk(i,1:length(k))=k';
end
figure(1)
% subplot(131)
VPest=log(abs(pp(:,1)))*Fx;
plot(fftshift(fff(1,1050:1175)),VPest(1050:1175,1),'r.',fftshift(fff(1,1325:1450)),VPest(1325:1450,1),'r.');hold off;
xlabel('Fr¨¦quence (Hz)')
ylabel('Pole')
title('Pole')
figure(2)

VPest=abs(2*pi*fff'./(fftshift(angle(pp(:,1))*Fx)+1e-3));
for rn = 1:2500
    if(VPest(rn,1)>1e5)
        VPest(rn,1)=0;
    end
end
%VPest=abs(2*pi*fff'./fftshift(angle(pp(:,1))*Fx));
semilogy(fff(1,1050:1175),VPest(1050:1175,1),'r.',fff(1,1325:1450),VPest(1325:1450,1),'r.');hold off;
xlabel('Fr¨¦quence (Hz)')
ylabel('Vitesse de phase (m/s)')
title('vvx')
%%

figure(3)
fftx=fftshift(fft(imageox,[],2));
imagesc(real(fftx));
xlabel('Fr¨¦quence (Hz)')
ylabel('Coordonne')
title('fft-x')


