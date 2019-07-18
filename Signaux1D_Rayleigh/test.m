
im= imread('test1.png');
im2=im2double(im);
figure(1)
subplot(2,3,1)
imshow(im2)
title('image')
subplot(2,3,2)
imshow(im(:,:,1))
title('image R')
subplot(2,3,3)
imshow(im(:,:,2))
title('image G')
subplot(2,3,4)
imshow(im(:,:,3))
title('image B')
%%
sigma = 2;
subplot(2,3,5)
RR = im(:,:,1);
thresh=[0.01,0.17];
f=edge(double(RR),'canny',thresh,sigma);
imshow(f)
title('Bord ')
[H,T,R]=hough(f,'ThetaResolution',45,'RhoResolution',10);
P=houghpeaks(H,400,'Threshold',80,'NHoodSize',[1,1]);

subplot(2,3,6)
lines = houghlines(f,T,R,P,'FillGap',3,'Minlength',8);
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    len = norm(lines(k).point1 - lines(k).point2);
    Len(k)=len;
    if(len>50 & len<200)
        plot(xy(:,1),-xy(:,2),'LineWidth',1,'color','r');
        hold on;
    end
end
title('Hough ')
%%
% figure(3)
% I=im;
% 
% Ihsv=rgb2hsv(I);
% Iv=Ihsv(:,:,1);                    %提取v空间
% Ivl=Iv(500:end,:);              %截取下半部
% Iedge=edge(Ivl,'sobel');    %边沿检测
% Iedge = imdilate(Iedge,ones(3));%图像膨胀
% 
% %新建窗口，绘图用
% imshow(Iedge);
% hold on
% 
% %左方直线检测与绘制
% %得到霍夫空间
% %[H1,T1,R1] = hough(Iedge,'Theta',20:0.1:75);
% [H1,T1,R1] = hough(Iedge,'Theta',-75:0.1:-20);
% 
% %求极值点
% Peaks=houghpeaks(H1,5);
% 
% %得到线段信息
% lines=houghlines(Iedge,T1,R1,Peaks);
% 
% %绘制线段
%  for k=1:length(lines)
% xy=[lines(k).point1;lines(k).point2];   
% plot(xy(:,1),xy(:,2),'LineWidth',4);
%  end
% 
%  
%  %右方直线检测与绘制
% [H2,T2,R2] = hough(Iedge,'Theta',-75:0.1:-20);
% Peaks1=houghpeaks(H2,5);
% lines1=houghlines(Iedge,T2,R2,Peaks1);
% for k=1:length(lines1)
% xy1=[lines1(k).point1;lines1(k).point2];   
% plot(xy1(:,1),xy1(:,2),'LineWidth',4);
% end
% 
% hold off
