%% EstDStates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This m-file provides diffusion coefficient estimation using 2-dimensional 
% trajectories of Brownian particles such as single fluorescent molecules 
% observed under TIRF microscopy. 
%
% The estimation is based on not mean squared displacement but statistic 
% distribution of displacement, which enables an estimation of individual 
% diffusion coefficients even if the trajectories are heterogeneous. 
%
% The most likely number of diffusion states with different diffusion 
% coefficients is suggested based on Akaike Information Criterion. 
% 
% It is optional to analyze transition kinetics among the diffusion states.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User data
% Please define the file name, scale and dt below. 

XY=csvread('data.csv'); % single molecule trajectory data
scale=0.08; % micrometer/pixel if needed 
dt=1/30; % time interval of measurement

%% Data formatting

T=0:dt:dt*max(XY(:,1)); % actual time (sec) along the trajectory

i=0;
for k=1:length(XY)
    frame=XY(k,1);
    if frame==0 
        R=zeros(length(T), 2);
        i=i+1;
        j=0;
    end
    j=j+1;
    temp(j,1)=XY(k,2)*scale;
    temp(j,2)=XY(k,3)*scale;
    A{1,i}=temp;
    BD(1,i)=j;
end    
N=i; % number of molecules
clear frame temp

%% Displacement calculation
% Displacement, Dr, that a molecule moved within a unit time interval, Dt,
% will be calculated. For the purpose of discriminating multiple diffusion 
% states, it is recommended to use the shortest Dt as possible, that is,
% Dt=fint*dt (fint=1).

fint=1;

for i=1:N
    F2=(length(T)-1)*(i-1)+1;F3=(length(T)-1)*i;
    SD=zeros(length(T)-1,length(T)-1);
    bd=BD(1,i);
    Nd=bd-1:-1:1;
    if bd==1
    else
        for j=1:bd-1 % lagtime (frame) for displacement measurement
            for k=1:bd-j % time (frame) along the trajectory
                SD(k,j)=(A{1,i}(k+j,1)-A{1,i}(k,1))^2+(A{1,i}(k+j,2)-A{1,i}(k,2))^2;
            end
            temp(:,1)=SD(:,j);
            MSD(j,i)=sum(temp)/Nd(1,j);
        end
        R=sqrt(SD);
        Dr(F2:F3,1)=R(:,fint);
        Dr_T(:,i)=R(:,fint); 
    end
    clear Nd
end
clear F2 F3 SD bd temp R

Dr(Dr==0,:)=[]; % dataset of Dr

%% Localization error estimation
% The position of a molecule in xy coordinates determined by fitting the 
% point spread function of the tagged fluorophore contains an error. 
% The variance, V, of the localization error follows MSD(t)=4*D*t+4*V, 
% where D is an averaged diffusion coefficient. MSD values calculated with 
% as short lagtimes as possible are used for the extrapolation to exclude 
% effects of super- or sub-diffusion.

y=3; % the trajectories longer than y frames are used for MSD calculation
z=2; % MSD(1*dt) to MSD(z*dt) are used for the extrapolation

MSDmean=zeros(1,y);
n=0;
for i=1:N
    if BD(1,i)>y
        MSDmean(1,:)=MSDmean(1,:)+MSD(1:y, i)';
        n=n+1;
    else
    end
end
MSDmean=MSDmean/n;
clear n

params0 = [0.01,0.001]; % initial guess for D and V
ub=[inf,inf];
lb=[0,0];
options=optimset('MaxFunEvals',1000,...
    'MaxIter',10000,...
    'TolFun',1e-8,...
    'TolX',1e-7);
[params,resnorm] = lsqcurvefit(@modelMSDFunction,params0,...
    T(1,2:z+1),MSDmean(1,1:z),lb,ub,options);

figure,plot(T(1,2:z+1),MSDmean(1,1:z),'r+');
hold on;
plot(T(1,1:z+1),modelMSDFunction(params,T(1,1:z+1)),'k:');
xlabel('Time (sec)');
ylabel('Mean square displacement');
legend('data','fitting');
hold off

Vest=params(2); % estimate of the variance

%% Diffusion coefficient estimation
% Estimate diffusion coefficients by maximum likelihood estimation (MLE). 
% MLE will be repeated 4 times assuming 4 different models in which the 
% number of molecular states with different diffusion coefficients is 1, 
% 2, 3 or 4. Sometimes better estimation will be provided by changing the 
% initial guess. Akaike Information Criterion (AIC) will be the minimum 
% when the model is the most likely model. 

Dr(1,2)=Vest;Dr(2,2)=dt*fint;Dr(3,2)=length(Dr);

% MLE assuming single state
params0=0.1; % initial guess for D
lb=0;
ub=inf;
options=optimset('MaxFunEvals',1000,...
        'MaxIter',10000,'TolFun',1e-15,'TolX',1e-10,'Display','iter');
[params,fval]=fmincon(@(params)model1stateFunction(params,Dr),...
        params0,[],[],[],[],lb,ub,[],options);
Dest(1,1)=params;Dest(5,1)=1;
L(1,1)=fval;

% MLE assuming two states
params0=[0.02;0.6;0.3]; % initial guess for D1, D2 and fraction of D1
lb=[0;0;0];
ub=[inf;inf;1];
options=optimset('MaxFunEvals',1000,...
        'MaxIter',10000,'TolFun',1e-15,'TolX',1e-10,'Display','iter');
[params,fval]=fmincon(@(params)model2stateFunction(params,Dr),...
        params0,[],[],[],[],lb,ub,[],options);
Dest(1:2,2)=params(1:2,1);
Dest(5,2)=params(3,1);
Dest(6,2)=1-params(3,1);
L(1,2)=fval;

% MLE assuming three states
params0=[0.01;0.05;0.5;0.3;0.5;0.2]; % initial guess for D1, D2, D3 and their fraction
lb=[0;0;0;0;0;0];
ub=[inf;inf;inf;inf;inf;inf];
options=optimset('MaxFunEvals',5000,...
        'MaxIter',10000,'TolFun',1e-15,'TolX',1e-10,'Display','iter');
[params,fval]=fmincon(@(params)model3stateFunction(params,Dr),...
        params0,[],[],[],[],lb,ub,[],options);
Dest(1:3,3)=params(1:3,1);
Dest(5,3)=params(4)/(params(4)+params(5)+params(6));
Dest(6,3)=params(5)/(params(4)+params(5)+params(6));
Dest(7,3)=params(6)/(params(4)+params(5)+params(6));
L(1,3)=fval;

% MLE assuming four states
params0=[0.01;0.05;0.10;0.60;0.001;0.001;0.001;0.001]; % initial guess
lb=[0;0;0;0;0;0;0;0];
ub=[inf;inf;inf;inf;inf;inf;inf;inf];
options=optimset('MaxFunEvals',5000,...
        'MaxIter',10000,'TolFun',1e-15,'TolX',1e-10,'Display','iter');
[params,fval]=fmincon(@(params)model4stateFunction(params,Dr),...
        params0,[],[],[],[],lb,ub,[],options);
Dest(1:4,4)=params(1:4,1);
Dest(5,4)=params(5)/(params(5)+params(6)+params(7)+params(8));
Dest(6,4)=params(6)/(params(5)+params(6)+params(7)+params(8));
Dest(7,4)=params(7)/(params(5)+params(6)+params(7)+params(8));
Dest(8,4)=params(8)/(params(5)+params(6)+params(7)+params(8));
L(1,4)=fval;

% to display a histogram of displacement
drbin=0.002; % size of bin 
drmax=1.0; % max value of displacement in the histogram
dr=drbin/2:drbin:drmax;

Drdist=hist(Dr(:,1),dr)/length(Dr)/drbin;
MLE(1,:)=1/(2*Dest(1,1)*Dr(2,2)+2*Vest)*dr.*exp(-dr.*dr/...
    (4*Dest(1,1)*Dr(2,2)+4*Vest));
MLE(2,:)=Dest(5,2)/(2*Dest(1,2)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(1,2)*Dr(2,2)+4*Vest))...
    +Dest(6,2)/(2*Dest(2,2)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(2,2)*Dr(2,2)+4*Vest));
MLE(3,:)=(Dest(5,3)/(2*Dest(1,3)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(1,3)*Dr(2,2)+4*Vest))...
    +Dest(6,3)/(2*Dest(2,3)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(2,3)*Dr(2,2)+4*Vest))...
    +Dest(7,3)/(2*Dest(3,3)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(3,3)*Dr(2,2)+4*Vest)));
MLE(4,:)=(Dest(5,4)/(2*Dest(1,4)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(1,4)*Dr(2,2)+4*Vest))...
    +Dest(6,4)/(2*Dest(2,4)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(2,4)*Dr(2,2)+4*Vest))...
    +Dest(7,4)/(2*Dest(3,4)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(3,4)*Dr(2,2)+4*Vest))...
    +Dest(8,4)/(2*Dest(4,4)*Dr(2,2)+2*Vest)*...
    dr.*exp(-dr.*dr/(4*Dest(4,4)*Dr(2,2)+4*Vest)));

figure,semilogx(dr,MLE(1,:),'r');
hold on;
plot(dr,MLE(2,:),'g');
plot(dr,MLE(3,:),'b');
plot(dr,MLE(4,:),'c');
plot(dr,Drdist,'o');
legend('1 state','2 states','3 states','4 states','data');
xlabel('Displacement (um)');
ylabel('Probability density');
hold off;

penalty=[1;3;5;7];
AIC(1,1)=2*L(1,1)*Dr(3,2)+penalty(1,1)*log(log(Dr(3,2)));
AIC(2,1)=2*L(1,2)*Dr(3,2)+penalty(2,1)*log(log(Dr(3,2)));
AIC(3,1)=2*L(1,3)*Dr(3,2)+penalty(3,1)*log(log(Dr(3,2)));
AIC(4,1)=2*L(1,4)*Dr(3,2)+penalty(4,1)*log(log(Dr(3,2)));

figure,plot(AIC,'b');
hold on;
xlabel('Number of states');
ylabel('AIC');
hold off;

%% State transition kinetics (optional)
% The number of molecules adopting each diffusion state as a function of 
% time after an onset of the trajectory will be calculated as well as 
% the time series of the fraction of the states. An example for the case 
% of 3 states is described below.

BDbin=1:1:length(T);
BDhist=hist(BD(1, :), BDbin);
temp(1, 1)=N;
for l=1:length(T)-1
    temp(1, l+1)=temp(1, l)-BDhist(1, l);
end
RC=temp/N;
clear temp

s=90; % maximum time (frame) of analyzing transition kinetics

for i=1:s
    temp=Dr_T(i,:)';temp(temp==0)=[];
    temp(1,2)=Vest;temp(2,2)=dt*fint;temp(3,2)=length(temp);
    temp(4,2)=Dest(1,3);temp(5,2)=Dest(2,3);temp(6,2)=Dest(3,3);
    
    params0 = [0.5;0.45;0.1];
    lb = [0;0;0];
    ub = [1;1;1];
    options=optimset('MaxFunEvals',1000,...
        'MaxIter',10000,'TolFun',1e-15,'TolX',1e-15,'Display','iter');
    [params,fval]=fmincon(@(params)model3state2Function(params,temp),...
        params0,[],[],[],[],lb,ub,[],options);
    F(1,i)=params(1)/(params(1)+params(2)+params(3));
    F(2,i)=params(2)/(params(1)+params(2)+params(3));
    F(3,i)=params(3)/(params(1)+params(2)+params(3));   
    RC(2:4,i+1)=F(1:3,i)*RC(1,i+1);
    clear temp
end

figure,plot(T(1,2:s+1),RC(2,2:s+1),'b+');
hold on
plot(T(1,2:s+1),RC(3,2:s+1),'g+');
plot(T(1,2:s+1),RC(4,2:s+1),'r+');
plot(T(1,1:s+1),RC(1,1:s+1),'k');
legend('D1','D2','D3','all');
xlabel('Time (sec)');
ylabel('Number of molecules');
hold off

figure,plot(T(1,2:s+1),F(1,1:s),'b+');
hold on
plot(T(1,2:s+1),F(2,1:s),'g+');
plot(T(1,2:s+1),F(3,1:s),'r+');
legend('D1','D2','D3');
xlabel('Time (sec)');
ylabel('Ratio');
hold off
