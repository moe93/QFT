function [ tfobj ] = sym2tf( symobj, Ts)
% SYM2TF convert symbolic math rationals to transfer function

if isnumeric(symobj)
    tfobj=symobj;
    return;
end

[n,d]=numden(symobj);
num=sym2poly(n);
den=sym2poly(d);

if nargin==1
    tfobj=tf(num,den);
else
    tfobj=tf(num,den,Ts);
end