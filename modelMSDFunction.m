function F = modelMSDFunction(parameters,xdata)

D=parameters(1);
a=parameters(2);

F=4*D*xdata+4*a;