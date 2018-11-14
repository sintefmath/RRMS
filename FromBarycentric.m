%%%Copyright 2018 SINTEF AS
function x = FromBarycentric(a,b,c,d, lambda)
    A = [[a;b;c;d]';[1,1,1,1]];
    x = A*lambda;
    x = x(1:3)';
end