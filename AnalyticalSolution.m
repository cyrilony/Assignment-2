function V = AnalyticalSolution(xInd,yInd,Vo,xMax,yMax,nMax)

V1 = (4*Vo)/pi;
a=yMax;
b=xMax/2;
x=xInd-b;
y=yInd;
nList = 2*(0:round((nMax-1)*0.5))+1;
V = V1*sum((sin(nList*(pi*y/a))./(nList)).*(...
    cosh(nList*(pi*(x)/yMax))./cosh(nList*(pi*b/a))));
end

