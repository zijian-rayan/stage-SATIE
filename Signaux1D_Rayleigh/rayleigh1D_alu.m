%Calcul des param�tres:
clear;clc;
load config_file
load signal_600_1.mat

tmp=resample(T,1,10);
nt=length(tmp)-1;
dt=median(diff(tmp));
t0=T(1);
t=t0+(0:nt-1)*dt;

nx=600;
ny=1;
dx=0.2e-3;
dt=median(diff(t));
x=(0:nx-1)*dx;
k=(((0:nx-1)/nx)-.5)/dx;
f=(((0:nt-1)/nt)-.5)/dt;


%%
%Chargement des fichiers dans la matrice s(x,y,t)
load s

%Affichage s:
figure(1)
subplot(2,3,1)
imagesc(t,x,s),colorbar,hold,plot(t,t*2900),hold
xlabel('Temps (s)')
ylabel('Position (m)')
title('Onde de Rayleigh dans l''axe du transducteur: Signaux s(x,t)')
%%
subplot(2,3,2)
%Spectre 2D et comparaison avec th�orie:
imagesc(f,-k,fftshift(abs(fft2(s)))),hold,plot(f,f./2930),hold
xlabel(' f ')
ylabel('k ')
title('Spectre 2D et comparaison')
%% dsp
subplot(2,3,3)
y=s;   % ��ֵ��λ Ϊ1nm��
LONG=30e-6;  %LONG����30e-6 m;                                                            Ϊ0.03m*0.03m  �����PERSSON   С��3��������
[L, L1] = size(y);    % ��Ϊ�� ������ ����L=L1   LΪ�ܵĸ���      ȡż����������
a=LONG/L; % ����ָ��С���    a= LONG/L;  �� L ��1024����                                  ԼΪ 3e-5 m�� ��PERSSON 1e-4      ��΢Сһ��㣬  ԼΪ1/3��
Fs=1/a;      %����Ƶ�ʣ�  Hz
% mesh(y)
% title('mesh(y)')


% subplot(2,3,4)
% surf(y)
% camlight left; lighting phong
% shading interp
% title('surf(y)')
MEAN=mean(mean(y));    
EE=1/(L^2)*sum(sum((y-MEAN).^2));    
FF=fft2(y);          
ABSFF=abs(FF);
dx = a;            
dy = a;
g=y;
g0 = fftshift(g);     % shift
G0 = fft2(g0)*dx*dy;  % 2D fft and dxdy scaling
G = fftshift(G0);     % center                          �൱�ڸ���Ҷ�仯
ABSFF2=abs(G)/(dx*dy);     % ABSFF1=abs(G)*Fs^2;
ABSFF2((L/2+1),(L/2+1))=1/4*(ABSFF2((L/2+2),(L/2+1))+ABSFF2((L/2),(L/2+1))+ABSFF2((L/2+1),(L/2+2))+ABSFF2((L/2+1),(L/2)));   %Ŀ�ļ�ȥ����㡣
fx = -1/(2*dx):1/(L*a):1/(2*dx)-(1/(L*a)); % x freq coords                   % Ƶ�ʷֲ�
fy = -1/(2*dy):1/(L*a):1/(2*dy)-(1/(L*a)); % y freq coords

% subplot(2,3,5)
[fx,fy] = meshgrid(fx,fy);
xx = zeros(600,600);
% for i = 1:600
%     for j = 101:700
%         jj=j-100;
%         xx(i,jj) = ABSFF2(i,j);
%     end
% end
% surf(fx,fy,xx)  % display transform magnitude
% surf(ABSFF2)
% camlight left; lighting phong
% shading interp
% ylabel('freq. fy (1/m)'); xlabel('freq.fx (1/m)');zlabel('dsp C (m^2*m^2) ')
% title('magnitude');
Pxx1=ABSFF.^2/(L^2*Fs^2);   %���persson�� ������һ��  |FF|^2*(a^2)/A;
Pxx2=ABSFF2.^2/(L^2*Fs^2);

subplot(2,3,6)

for i = 1:600
    for j = 101:700
        jj=j-100;
        xx(i,jj) = Pxx2(i,j);
    end
end
surf(fx,fy,xx)  % display transform magnitude
surf(Pxx2)
camlight left; lighting phong
shading interp
ylabel('freq. fy (1/m)'); xlabel('freq.fx (1/m)');zlabel('dsp C (m^2*m^2) ')
title('magnitude');
%%  autocorelation
figure(1)
xco=xcorr2(s);
dsp = fftshift(abs(fft2(xco)));
subplot(2,2,1)
imagesc(xco);
title('image Auto-corr��lation');
subplot(2,2,2)
imagesc(f,-k,fftshift(abs(fft2(xco)))),hold,plot(f,f./2930),hold
title('image dsp');
subplot(2,2,3)
surf(xco)
title('surf Auto-corr��lation');
camlight left; lighting phong
shading interp
subplot(2,2,4)
surf(fftshift(abs(fft2(xco))))
title('surf dsp');
camlight left; lighting phong
shading interp
%% dsp -> ifft
figure(3)
subplot(2,2,2)
imagesc(Pxx2);
title('image dsp');
subplot(2,2,1)
imagesc(f,-k,fftshift(abs(ifft2(Pxx2)))),hold,plot(f,f./2930),hold
title('image Auto-corr��lation');
subplot(2,2,4)
surf(Pxx2)
title('surf dsp');
camlight left; lighting phong
shading interp
subplot(2,2,3)
surf(fftshift(abs(ifft2(Pxx2))))
camlight left; lighting phong
shading interp
title('surf Auto-corr��lation');

%%

c=zeros(1199,2001);
for i = 550 : 650
    for j = 950 : 1050
        c(i,j)=1;
    end
end
for i = 598 : 602
    for j = 999 : 1003
        c(i,j)=0;
    end
end
subplot(2,2,1)
surf(xco)
title('surf Auto-corr��lation');
camlight left; lighting phong
shading interp
subplot(2,2,2)
ff = fftshift(abs(fft2(xco)));
surf(ff)
title('surf dsp');
camlight left; lighting phong
shading interp
subplot(2,2,3)
ff=ff .* c;
surf(ff)
title('surf dsp');
camlight left; lighting phong
shading interp
subplot(2,2,4)
surf(fftshift(abs(ifft2(ff))))
title('surf dsp');
camlight left; lighting phong
shading interp







