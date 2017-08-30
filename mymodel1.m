function [newS, newU, newB]=mymodel1(mygrid,params,pp1,pp2,pp3,L);

		

params_S=params(1:pp1);
params_B=params((pp1+1):(pp1+pp2));
params_U=params((pp1+pp2+1):(pp1+pp2+pp3));



% Basal model
C=zeros(1,sum(L)-1024);
C(1:32)=params_B;
newB=waverec(C,L,'db3');




% surface model
a=params_S(1);
b=params_S(2);
c=params_S(3);

newS=a*mygrid.^2 + b*mygrid + c;


% velocity model
rho=911;
g=9.81;

%dx=mygrid(2)-mygrid(1);

%dSS = diff(newS);
%dSS=[dSS,dSS(length(dSS))];
%dSS=(1/dx)*dSS;

dSS = 2*a*mygrid + b;

HH=newS-newB;

u_b = params_U(1);
a0=params_U(2);
a1=params_U(3);
b1=params_U(4);
a2=params_U(5);
b2=params_U(6);
a3=params_U(7);
b3=params_U(8);
w=params_U(9);

x=mygrid;

A=a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w);

newU = u_b + A.*((0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3);
