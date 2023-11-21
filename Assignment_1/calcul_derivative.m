clear all;clc;

syms x y z mu ms t rho omegas mus;

omega4=0.5*(x^2+y^2)+(1-mu)/sqrt((x+mu)^2+y^2)+...
mu/sqrt((x+mu-1)^2+y^2)+0.5*mu*(mu-1)+ms/sqrt((x-rho*cos(omegas*t))^2+ ...
(y-rho*sin(omegas*t))^2)-(ms/rho^2)*(x*cos(omegas*t)+y*sin(omegas*t));

omega4x=diff(omega4,x)
omega4y=diff(omega4,y)

omega4xx=diff(omega4x,x);
omega4xy=diff(omega4x,y);
omega4yy=diff(omega4y,y);
omega4yx=diff(omega4y,x);