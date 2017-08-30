clear all;
load running_diff_langevin;

s0=10;
e0=sample_count;

figure(1);
subplot(311)
plot(paramskeep(s0:e0,1));
subplot(312)
plot(paramskeep(s0:e0,15));
subplot(313)
plot(paramskeep(s0:e0,37));




if 1
	figure(3);
	plot(paramskeep(s0:e0,36)); title('u_b');

	figure(4);
	%plot(paramskeep(s0:e0,37));
	for i=1:8
		subplot(4,2,i);
		plot(paramskeep(s0:e0,36+i));
	end
end




if 0
	s0=1;
	e0=sample_count;

	idx=s0:e0;
	mp=mean(paramskeep(idx,:));



	[meanS, meanU, meanB] = mymodel1(mygrid,mp,pp1,pp2,pp3,L); 


	figure(10);
	plot(mygrid,SS,'r');hold on;
	plot(mygrid, meanS,'k'); hold off;

	figure(11);
	plot(mygrid, UU, 'r');hold on;
	plot(mygrid, meanU,'k');hold off;



	figure(12);
	newB=zeros(1,1024);

	C=zeros(1,sum(L)-1024);
	C(1:32)=mp((pp1+1):(pp1+pp2));
	meanB=waverec(C,L,'db3');
	plot(BB);hold on;
	plot(meanB,'r');hold off;
end