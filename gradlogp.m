function glogp = gradlogp(mygrid, SS,BB,UU, params, mu_prior, theta_prior, L, W, Jd, pp1, pp2, pp3);




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

newS=a*mygrid.^2+ b*mygrid+ c;


% velocity model
rho=911;
g=9.81;

dSS = 2*a*mygrid+b;

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


		
p=length(params);
%gradprior=zeros(1,p);
gradprior= - (params-mu_prior)./theta_prior;






n=1024;



gradS=zeros(1,p);
GS=1;
if GS
	%===================================================================================================
	% surface model

	dSda = mygrid.^2; 
	dSdc = ones(1,n);
	dSdb = mygrid;


	gradS(1) = sum(dSda.*(newS-SS));
	gradS(2) = sum(dSdb.*(newS-SS));
	gradS(3) = sum(dSdc.*(newS-SS));

	gradS = gradS * (-2) * Jd(1);
	%===================================================================================================
end



gradB=zeros(1,p);
GB=1;
if GB	
	gradB((pp1+1):(pp1+pp2)) = (newB-BB)*W';

	gradB = gradB * (-2) *Jd(2);
	%===================================================================================================
end


gradU=zeros(1,p);
GU=1;
if GU
	%===================================================================================================
	% surface model

	dSSda = 2*mygrid;
	T1 = 4*dSda.*(HH.^3).*(dSS.^3);
	T2 = (HH.^4) .* 3 .* dSSda.*(dSS.^2);
	dUda = (0.5)*A*(rho * g)^3.*(T1+T2);

	dSSdb = ones(1,n);
	T1 = 4*dSdb.*(HH.^3).*(dSS.^3);
	T2 = (HH.^4) .* 3.* dSSdb .* (dSS.^2);
	dUdb = (0.5)*A*(rho * g)^3.*(T1+T2);

	dSSdc = zeros(1,n);
	T1 = 4*dSdc.*(HH.^3).*(dSS.^3);
	T2 = zeros(1,n);
	dUdc = (0.5)*A*(rho * g)^3.*(T1+T2);



	gradU(1) = sum(dUda.*(newU-UU));
	gradU(2) = sum(dSdb.*(newU-UU));
	gradU(3) = sum(dSdc.*(newU-UU));



%	gradU((pp1+1):(pp1+pp2)) = -(  0.5*A*(rho*g)^3*4.*(HH.^3).*(newU-UU) ) * W';
%	for j=1:pp2
%		dUdcj = - (4/2)*(rho*g)^3.*(HH.^3).*(dSS).^3.*W(j,:);
%		gradU(pp1+j) = sum( dUdcj.* (newU-UU) ); 
%	end


	dUdub = ones(1,n);
	gradU(pp1+pp2+1) = sum(dUdub.*(newU-UU));


	dUda0 = (0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	dUda1 = cos(x*w).*(0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	dUdb1 = sin(x*w).*(0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	dUda2 = cos(2*x*w).*(0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	dUdb2 = sin(2*x*w).*(0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	dUda3 = cos(3*x*w).*(0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	dUdb3 = sin(3*x*w).*(0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	
	dUdA = (0.5)*(rho*g)^3.*(HH.^4).*(dSS).^3;
	dAdw = -a1.*x.*sin(x*w) + b1.*x.*cos(x*w) - 2*a2.*x.*sin(2*x*w) + 2*b2.*x.*cos(2*x*w) - 3*a3.*x.*sin(3*x*w) + 3*b3.*x.*cos(3*x*w);
	
	
	dUdw = dUdA.*dAdw;
	
	
	gradU(pp1+pp2+2) = sum(dUda0.*(newU-UU));
	gradU(pp1+pp2+3) = sum(dUda1.*(newU-UU));
	gradU(pp1+pp2+4) = sum(dUdb1.*(newU-UU));
	gradU(pp1+pp2+5) = sum(dUda2.*(newU-UU));
	gradU(pp1+pp2+6) = sum(dUdb2.*(newU-UU));
	gradU(pp1+pp2+7) = sum(dUda3.*(newU-UU));
	gradU(pp1+pp2+8) = sum(dUdb3.*(newU-UU));
	gradU(pp1+pp2+9) = sum(dUdw.*(newU-UU));

	



	gradU = gradU * (-2) * Jd(3);
	%===================================================================================================
end





glogp = gradprior + gradS + gradB + gradU;
%glogp(37)
%gradprior(36)
%gradS(36)
%gradB(36)
%gradU(36)
%glogp(7)