%%%%%%%%%%%% Rayleigh Analyse %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ARFAOUI MERZOUK %%%%%%%%%%%%%%%%%%%%%
clear all
% cd('C:\Users\SATIE\Desktop\Salim_Said\201120104_test_3_XY')
% cd('C:\DATA\work\recherche\stagiaires\2014_said_arfaoui_salim_merzouk_rayleigh\Salim_Said\201120104_test_3_XY')
load config_file
load signal_1_1.mat
tmp=resample(T,1,10);
nt=length(tmp);
t0=T(1);
dx=0.2e-3
dy=0.2e-3
dt=median(diff(tmp));

%Nombre de Point
nx=1; %sur l'axe X
ny=1; % sur l'axe Y

%Chargement des fichiers dans la matrice s(x,y,t)
s=zeros(nx,ny,nt);
for j=1:ny,
    for i=1:nx,
    eval(sprintf('load(''signal_%d_%d.mat'')',i,j));
    tmp=resample(Y,1,10);
    s(i,j,:)=tmp;
    i,j
    end
end
s=squeeze(s);
save s s

%Calcul des paramtètres et affichage de la matrice s(x,y,t)
[nx,ny,nt]=size(s);
x=(0:nx-1)*dx;
y=(0:ny-1)*dy;
t=t0+(0:nt-1)*dt;
%t=tmp;
k=(((0:nx-1)/nx)-.5)/dx;
f=(((0:nt-1)/nt)-.5)/dt;
% afichage de la courbe s(x,y,t)
figure (1)
imagesc(t,x,squeeze(s)),colorbar,hold
plot(t,t*3177,'r'),hold
xlabel('Temps (s)')
ylabel('Position (m)')
title('suivi de l''onde de Rayleigh sur XY: signaux s(x,y,t) Y = 50')

%Spectre 2D
figure(2)
imagesc(f,-k,log(fftshift(abs(fft2(squeeze(s(50)))))))%,hold,plot(f,f./2930), hold
xlabel('Frequence (Hz)')
ylabel('Vecteur d''onde (m^-1)')
title('Onde de Rayleigh sur XY: courbes de dispersion, pour Y=50')

vidObj = VideoWriter('Rayleigh_XY.avi');
open(vidObj);



j=1,clear MM, mmm=max(s(:));
for i=100:200
    imagesc(squeeze(s(i)),[-mmm mmm])
    axis equal tight,colorbar 
    MM(j) = getframe
    j=j+1
    i
currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);
   end
close(vidObj);

% mmm=max(s(:));for i=100:1000,surf(squeeze(s(:,:,i))),axis equal tight,colorbar,i,pause,end

