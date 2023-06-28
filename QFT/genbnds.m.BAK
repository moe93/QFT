function bds = genbnds(ptype,w,W,A,B,C,D,Pnom,ph_d)
% GENBNDS Compute General QFT bounds.
%         GENBNDS(PTYPE,W,Ws,A,B,C,D,Pnom,PHS) computes bounds for the
%         configuration designated by PTYPE.  W is the set of frequencies
%         designating bounds to compute.  Ws is the performance specification
%         (scalar or LTI/FRD model). A, B, C and D are scalers and/or LTI/FRD arrays,
%         Pnom is the nominal plant needed to assign bounds to a nominal loop.
%         PPHS specifies the phase (degrees) grid for bound computations
%         (0 to -360 every -5 deg is default).
%
%         PTYPE specifies the closed-loop configurations.
%         PTYPE=10: |(A+BG)/(C+DG)|<Ws
%         PTYPE=11: (|A|+|BG|)/|(C+DG)|<Ws
%
%         See also SISOBNDS, GRPBNDS, PLOTBNDS, SECTBNDS.

% Author: Craig Borghesani
% Date: 5/21/94
% Revised: 2/16/96 1:06 PM V1.1 updates, 7/3/03 9:58AM v2.5 updates.
% Copyright (c) 2002, Terasoft, Inc.

if ptype < 10 | ptype > 11,
 error('GENBNDS only accepts problem types between 10 and 11');
end

if isa(A,'lti'),
   if length(w) == 1 & prod(size(A)) > 1,
      A = squeeze(freqresp(A, w));
   else
      A = squeeze(freqresp(A, w)).';
   end;%if
   if any(any(isnan(A))),
      error('Frequency vector for A inconsistent with w.');
   end
end;%if isa(A,'lti')

if isa(B,'lti'),
   if length(w) == 1 & prod(size(B)) > 1,
      B = squeeze(freqresp(B, w));
   else
      B = squeeze(freqresp(B, w)).';
   end;%if
   if any(any(isnan(B))),
      error('Frequency vector for B inconsistent with w.');
   end
end;%if isa(B,'lti')

if isa(C,'lti'),
   if length(w) == 1 & prod(size(C)) > 1,
      C = squeeze(freqresp(C, w));
   else
      C = squeeze(freqresp(C, w)).';
   end;%if
   if any(any(isnan(C))),
      error('Frequency vector for C inconsistent with w.');
   end
end;%if isa(C,'lti')

if isa(D,'lti'),
   if length(w) == 1 & prod(size(D)) > 1,
      D = squeeze(freqresp(D, w));
   else
      D = squeeze(freqresp(D, w)).';
   end;%if
   if any(any(isnan(D))),
      error('Frequency vector for D inconsistent with w.');
   end
end;%if isa(D,'lti')

if isa(Pnom,'lti'),
   if length(w) == 1 & prod(size(Pnom)) > 1,
      Pnom = squeeze(freqresp(Pnom, w));
   else
      Pnom = squeeze(freqresp(Pnom, w)).';
   end;%if
   if any(any(isnan(Pnom))),
      error('Frequency vector for Pnom inconsistent with w.');
   end
end;%if isa(Pnom,'lti')

[rma,cma]=size(A);[rmc,cmc]=size(C);
[rmb,cmb]=size(B);[rmd,cmd]=size(D);
row = [rma,rmb,rmc,rmd];
col = [cma,cmb,cmc,cmd];
maxr = max(row); maxc = max(col);
test_P = ones(maxr,maxc);

% removal of wbd
wbd = [];

if nargin==8,
 [w,wbd,W,jk,jk,jk,jk,jk,jk,jk,ph_r,info]=bndsdef(w,wbd,W,test_P,[],[],[],[],[],ptype);
elseif nargin==9,
 [w,wbd,W,jk,jk,jk,jk,jk,jk,jk,ph_r,info]=bndsdef(w,wbd,W,test_P,[],[],[],[],ph_d,ptype);
else
 error('Improper number of inputs');
end

[rmw,cmw]=size(W);
maxr = max(rmw,maxr);
nbd=length(w);

if length(Pnom)~=nbd,
 error('Nominal plant information inconsistent with frequency vector');
end

if any(row~=1 & row~=maxr),
 error('Cannot have different numbers of rows in A, B, C & D');
end

if cma~=1 & cma~=nbd,
 error('Matrix A inconsistent with frequency vector');
end

if cmb~=1 & cmb~=nbd,
 error('Matrix B inconsistent with frequency vector');
end

if cmc~=1 & cmc~=nbd,
 error('Matrix C inconsistent with frequency vector');
end

if cmd~=1 & cmd~=nbd,
 error('Matrix D inconsistent with frequency vector');
end

if nargout == 1,
 eval(['bds=genbnd',int2str(ptype),'(w,wbd,W,A,B,C,D,Pnom,ph_r,info);'],'qfterror(1,info)');
else
 eval(['genbnd',int2str(ptype),'(w,wbd,W,A,B,C,D,Pnom,ph_r,info);'],'qfterror(1,info)');
end
