function [ze]=elliptic_plot(a,b,x0,y0,cang)
%[a    b     x0     y0     alpha]
%[.1   .1   -.5    -.5   cang   %

n=100;
z0=x0+sqrt(-1)*y0;
theta=(0:n-1)/(n-1)*2*pi;
z=exp(sqrt(-1)*theta);
x=real(z);
y=imag(z);
ze=z0+(x*a+sqrt(-1)*y*b)*exp(sqrt(-1)*cang/180*pi);
