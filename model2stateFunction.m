function F = model2stateFunction(parameters,Q)

D1=parameters(1);
D2=parameters(2);
a=parameters(3);
sigma2=Q(1,2);
t=Q(2,2);
N=Q(3,2);
Q=Q(:,1);

%P=a/(2*D1*t)*Q.*exp(-Q.*Q/(4*D1*t))...
%    +(1-a)/(2*D2*t)*Q.*exp(-Q.*Q/(4*D2*t));
P=a/(2*D1*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D1*t+4*sigma2))...
    +(1-a)/(2*D2*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D2*t+4*sigma2));

L=log(P);
F=-(sum(L)/N);