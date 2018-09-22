function F = model4stateFunction(parameters,Q)

D1=parameters(1);
D2=parameters(2);
D3=parameters(3);
D4=parameters(4);
k=parameters(5);
l=parameters(6);
m=parameters(7);
n=parameters(8);
sigma2=Q(1,2);
t=Q(2,2);
N=Q(3,2);
Q=Q(:,1);

a=k/(k+l+m+n);
b=l/(k+l+m+n);
c=m/(k+l+m+n);
d=n/(k+l+m+n);

%P=a/(2*D1*t)*Q.*exp(-Q.*Q/(4*D1*t))...
%    +b/(2*D2*t)*Q.*exp(-Q.*Q/(4*D2*t))...
%    +c/(2*D3*t)*Q.*exp(-Q.*Q/(4*D3*t))...
%    +d/(2*D4*t)*Q.*exp(-Q.*Q/(4*D4*t));
P=a/(2*D1*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D1*t+4*sigma2))...
    +b/(2*D2*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D2*t+4*sigma2))...
    +c/(2*D3*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D3*t+4*sigma2))...
    +d/(2*D4*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D4*t+4*sigma2));

L=log(P);
F=-(sum(L)/N);