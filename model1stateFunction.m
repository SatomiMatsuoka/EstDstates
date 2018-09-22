function F = model1stateFunction(parameters,Q)

D=parameters(1);
sigma2=Q(1,2);
t=Q(2,2);
N=Q(3,2);
Q=Q(:,1);

%P=1/(2*D*t)*Q.*exp(-Q.*Q/(4*D*t));
P=1/(2*D*t+2*sigma2)*Q.*exp(-Q.*Q/(4*D*t+4*sigma2));

L=log(P);
F=-(sum(L)/N);