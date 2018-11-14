%%%Copyright 2018 SINTEF AS
function lambda = ToBarycentric(a,b,c,d, x)
    A = [[a;b;c;d]';[1,1,1,1]];
    x = [x';1];
    lambda = A\x;
end