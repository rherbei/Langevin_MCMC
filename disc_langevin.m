clear all;
rand('state',0); %to reproduce results
randn('state', 0);
%format bank; % two digit display 
warning off;


% prior mean values
a=-2.583e-9;
b=0.004151;
c=1551;

%u_b=36;
u_b=41;

a0 =  -1.214e-17 ; %(-1.832e-16, 1.589e-16)
a1 =  -2.837e-16 ; %(-5.796e-16, 1.216e-17)
b1 =  -8.451e-17 ; %(-1.331e-16, -3.588e-17)
a2 =  -1.151e-16 ; %(-3.449e-16, 1.146e-16)
b2 =  -1.581e-16 ; %(-2.272e-16, -8.903e-17)
a3 =  -1.378e-18 ; %(-1.46e-16, 1.432e-16)
b3 =  -1.142e-16 ; %(-1.73e-16, -5.54e-17)
w  =   1.571e-05  ; %(1.392e-05, 1.75e-05)


% variance values
thetaS = 200;
thetaB = 2000;
thetaU = 200;


% load data set
load ../Data/icemodel1;
BB=YY;

% load wavelet coeffs
load ../Data/wmatrix;
W=W';

% surface parameters
params_S=[a,...		% -- 1
		  b,...		% -- 2
		  c];		% -- 3 

% basal wavelet coefficients
params_B=C(1:32);

% surface velocity parameters
params_U = [u_b,a0,a1,b1,a2,b2,a3,b3,w];



pp1=length(params_S);
pp2=length(params_B);
pp3=length(params_U);
pp=pp1+pp2+pp3;


params(1:pp1)=params_S;
params((pp1+1):(pp1+pp2))=params_B;
params((pp1+pp2+1):(pp1+pp2+pp3))=params_U;






		
MCS = 5000000;
MCO = 1000;


sample_count=1;
paramskeep=zeros(MCS/MCO+1,pp1+pp2+pp3);
thetakeep=zeros(MCS/MCO+1,3);

paramskeep(sample_count,1:pp1)=params_S;
paramskeep(sample_count,(pp1+1):(pp1+pp2))=params_B;
paramskeep(sample_count,(pp1+pp2+1):(pp1+pp2+pp3))=params_U;




	
%	priors means and variances
mu_prior=params;
theta_prior=(0.05*abs(params));



Jd=[1/(2*thetaS),...	% for surface data
	1/(2*thetaB),...	% for velocity data
	1/(2*thetaU)];	% for basal data




%step size settings 
h=zeros(1,pp);
hkeep=zeros(MCS/MCO+1,pp);

h(1)=1e-24;
%h(2)=1e-16;
%h(3)=1e-5;

h(4) = 1e-2;
h(5) = 1e-2;

h((pp1+1):(pp1+pp2)) = 1e-2*sqrt(theta_prior((pp1+1):(pp1+pp2)));

%-------------
% u_b
h(36)=1e-3;
%-------------
%a_0
h(37) = 1e-38;
%-------------
% a_1, b_1, a_2, b_2, a_3, b_3
h(38) = 1e-38;
h(39) = 1e-38;
h(40) = 1e-38;
h(41) = 1e-38;
%h(42) = 1e-39;
%h(43) = 1e-38;
%--------------
% w
%h(44) = 1e-19;
%%%%%%%%%%%%%%%%


% save grad-log-p
gkeep=zeros(MCS/MCO+1,pp);




%sample size
N=1024;

h_S=1e-5;
h_B=1e-3;
h_U=1e-5;


%gstart=gradlogp(mygrid, SS,BB,UU, params, mu_prior, theta_prior, L, W, Jd, pp1, pp2, pp3);


psiS=log(thetaS);
psiB=log(thetaB);
psiU=log(thetaU);

tic;

for k=1:MCS

	glogp=gradlogp(mygrid, SS,BB,UU, params, mu_prior, theta_prior, L, W, Jd, pp1, pp2, pp3);
	

	params = params + (h/2).*glogp + sqrt(h).*randn(1,pp);
	


	[newS, newU, newB]=mymodel1(mygrid,params,pp1,pp2,pp3,L);
	SSE_S=sum( (newS-SS).^2 );
	SSE_B=sum( (newB-BB).^2 );
	SSE_U=sum( (newU-UU).^2 );
	
	
	
	
	glog_psiS= exp(psiS)*(-(N/2)*(1/thetaS) + (SSE_S/2)*(1/thetaS^2));
	glog_psiB = exp(psiB)*(-(N/2)*(1/thetaB) + (SSE_B/2)*(1/thetaB^2));
	glog_psiU = exp(psiU)*(-(N/2)*(1/thetaU) + (SSE_U/2)*(1/thetaU^2));

	psiS=psiS+(h_S/2)*glog_psiS+sqrt(h_S)*randn;
	psiB=psiB+(h_B/2)*glog_psiB+sqrt(h_B)*randn;
	psiU=psiU+(h_U/2)*glog_psiU+sqrt(h_U)*randn;

	thetaS=exp(psiS);
	thetaB=exp(psiB);
	thetaU=exp(psiU);

	Jd=[1/(2*thetaS),...	% for surface data
		1/(2*thetaB),...	% for velocity data
		1/(2*thetaU)];		% for basal data


	
	if ~mod(k,MCO)
	
		

	 	
		disp(sprintf('%g  ', k));
		disp(sprintf('%g  ', params(36:40) ));
		disp(sprintf('%g  ', h(1:3) ));
		   
		disp(sprintf('%g  ', [thetaS thetaB thetaU]));
		%disp(sprintf('%g  ', [SSE_S/2, SSE_B/2 SSE_U/2]));	
      	disp(' ');
      	

	
		sample_count=sample_count+1;
		paramskeep(sample_count,:)=params;
		gkeep(sample_count,:)=glogp;
		thetakeep(sample_count,:)=[thetaS,thetaB,thetaU];

		if ~mod(k,50000)
			save running_diff_langevin paramskeep sample_count mygrid SS UU BB pp1 pp2 pp3 L h theta_prior Jd gkeep thetakeep;
		end
		
		toc
		tic

	end

	


end


save running_diff_langevin paramskeep sample_count mygrid SS UU BB pp1 pp2 pp3 L h theta_prior Jd gkeep thetakeep;
summarize_diff