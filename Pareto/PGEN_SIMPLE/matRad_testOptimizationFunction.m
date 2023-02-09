function f = matRad_testOptimizationFunction(x,maxPoints,pens)
    x= (x-maxPoints(1,:))./(maxPoints(2,:)-maxPoints(1,:));
    f = x*pens'; %x1*w1+x2*w2+x3*w3
    