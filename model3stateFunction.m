function F = model3stateFunction(parameters,Q)

D1=parameters(1);
D2=parameters(2);
D3=parameters(3);

k=parameters(4);
l=parameters(5);
m=parameters(6);
sigma2=Q(1,2);
t=Q(2,2);
N=Q(3,2);
Q=Q(:,1);

a=k/(k+l+m);
b=l/(k+l+m);
c=m/(k+l+m);

%P=a/(2*D1*t)*Q.*exp(-Q.*Q/(4*D1*t))...
%    +b/(2*D2*t)*Q.*exp(-Q.*Q/(4*D2*t))...
%    +c/(2*D3*t)*Q.*exp(-Q.*Q/(4*D3*t));

P=a/(2*D1*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D1*t+4*sigma2))...
    +b/(2*D2*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D2*t+4*sigma2))...
    +c/(2*D3*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D3*t+4*sigma2));

L=log(P);
F=-(sum(L)/N);