function [ symobj ] = tf2sym( tfobj, fp )
% TF2SYM convert transfer functions to symbolic math rationals

if isnumeric( tfobj )
    symobj = tfobj;
    return;
end

syms s;
[n, d] = tfdata( tfobj, 'v' );
num = poly2sym( n, s );
den = poly2sym( d, s );

if( nargin == 1 )
    symobj = num / den;
else
    symobj = vpa( num / den, fp );
end