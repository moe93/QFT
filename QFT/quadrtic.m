function [g1,g2] = quadrtic(a,b,c)
% QUADRTIC Solution of the quadratic. (Utility Function)
%          QUADRTIC computes the roots of the quadratic equations set up by
%          SISOBNDS and GENBNDS functions.  It also determines if there is
%          no possible linear, time-invariant controller that can solve the
%          problem.

% Author: Craig Borghesani
% Date: 1/22/93
% Revised: 2/16/96 10:45 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


[ra, ca] = size(a);
myeps    = 1e-16;

% set up matrices for states of variables.  It is these matrices that
% are used to determine whether there is no LTI controller that can
% solve the problem
stata = ( a <= 0 ); statb = ( b <= 0 ); statc = ( c <= 0 );
stat1 = ( a < 0  ); stat2 = ( b >  0 ); stat3 = ( c <  0 ); stat4 = (b.^2 < 4*a.*c);

lnega = find( a < 0 ); lzera = find( a == 0 ); lposa = find( a > 0 );
lnegb = find( b < 0 ); lposb = find( b >  0 );

cas2    = stata + statb + statc;
cas2b   = stat1 + stat2 + stat3 + stat4;

g1  = -ones( ra, ca );
g2  = g1;
g   = g1;

% --- If a == 0, then
%       g*b + c >= 0 ==> g >= -c/b
if( ~isempty(lzera) )
    locb1=[lposb(:); lnegb(:)];
    g( locb1 )      =-c( locb1 ) ./ b( locb1 );
    g1( lnegb )   = g( lnegb );
    g2( lposb )   = g( lposb );
end

loca1 = [ lnega(:); lposa(:) ];
g1( loca1 ) = (-b(loca1) + sqrt(b(loca1).^2-4*a(loca1).*c(loca1))) ./ (2*a(loca1));
g2( loca1 ) = (-b(loca1) - sqrt(b(loca1).^2-4*a(loca1).*c(loca1))) ./ (2*a(loca1));

z=find((g1<=0) | (imag(g1)~=0)); lz=length(z);
if lz==1,
    g1(z)=myeps;
elseif lz>1,
    g1(z)=ones(lz,1)*myeps;
end

z=find((g2<=0) | (imag(g2)~=0)); lz=length(z);
if lz==1,
    g2(z)=1/myeps;
elseif lz>1,
    g2(z)=ones(lz,1)/myeps;
end

k = (cas2==3 | cas2b==4);   % find all the bad guys

if ra > 1,
    bad_phases = any(k);  % determine only bad guy phases
else
    bad_phases = k;
end

% this row vector now tells us what individual phases are bad, so we
% replicate it
bad_phases = bad_phases(ones(ra,1),:);

% now replace the entire phase

%%%%%% V4.2 code
%g1(bad_phases) = g1(bad_phases)*0 + 248;
%g2(bad_phases) = g2(bad_phases)*0 - 248;

%%%%%% V5 code
% Reason: V5's boolean indexing change


loc_bad_phases = find(bad_phases);
g1(loc_bad_phases) = g1(loc_bad_phases)*0 + 248;
g2(loc_bad_phases) = g2(loc_bad_phases)*0 - 248;
