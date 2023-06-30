function bdb = sisobnd7(w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info)
%     Compute QFT bounds for the following closed-loop configuration
%
%                       |      PG     |
%                WS1 <= |F -----------| <= WS2
%                       |    1 + PGH  |
%
%     SISOBNDS(7,W,WBD,WS,P,R,NOM,C,LOC,PHS) computes bounds at frequencies designated by
%     WBD.  WS is the performance specification, P is the frequency response
%     data of the plant (complex), NOM designates the nominal plant and
%     controller, C. PHS specifies at which phases (degrees) to compute
%     bounds.
%
%     Note: R and LOC are fixed at 0 and 1, respectively
%
%     SISOBNDS(7,W,WBD,WS,P,[],C) computes bounds using default values
%     for NOM ([1,1]) and PHS ([0:-5:-360]).
%
%     WS for SISOBNDS(7,...) is different than all the other SISOBND# functions.
%     It is a matrix where WS1 represents the first row (upper spec) and WS2
%     represents the second row (lower spec).

% Authors: Yossi Chait and Craig Borghesani
% 6/17/92
% Copyright (c) 2003, Terasoft, Inc.


m=length(ph_r); nbd=length(w);
pbds=qsubset(wbd,w);
final_state=ones(1,nbd);
state2=[];
myeps=1e-16;

ones1=ones(1,m);

% preallocating matrix for more efficient programming
bmag=[ones(m,length(pbds))*myeps;ones(m,length(pbds))*1/myeps];

[rmp,cmp]=size(uP);[rpp,cpp]=size(vP);
[rmc,cmc]=size(uC);[rpc,cpc]=size(vC);
if nom(1)>rmp, error('Non-existent index for plant cases');
elseif nom(2)>rmc, error('Non-existent index for controller cases'); end

if rmc>1, error('Cannot have an uncertain controller.');
elseif rmp==1, error('Only one case in plant. No bounds are computed'); end

pct=1; cct=1; j=1;
while j<=length(pbds)

    if cmc>1, cct=pbds(j); end
    if cmp>1, pct=pbds(j); end

    phi = ph_r - (vP(nom(1),pct) + vC(nom(2),cct));
    for k=1:rmp-1
        mag_pk  = ones(rmp-k,1)*uP(k,pct);
        mag_pi  = uP(k+1:rmp,pct);
        ph_pk   = ones(rmp-k,1)*vP(k,pct);
        ph_pi   = vP(k+1:rmp,pct);
        phi     = phi(ones(rmp-k,1),:);
        psi_k   = ph_pk(:,ones1) + vC(cct) + phi;
        psi_i   = ph_pi(:,ones1) + vC(cct) + phi;
        cos_k   = cos(psi_k); cos_i=cos(psi_i);
        AA      = mag_pk.^2 .*mag_pi.^2 .*uC(cct)^2;
        Bk      = mag_pk.*mag_pi.*uC(cct);
        Bi      = Bk(:,ones1);
        delta   = W(pct)^2;
        A_k     = AA(:,ones1).*(1-1./delta);
        B_k     = Bk(:,ones1).*( mag_pk(:,ones1).*cos_i-mag_pi(:,ones1).*cos_k./delta);
        C_k     = mag_pk(:,ones1).^2-mag_pi(:,ones1).^2 ./delta;
        B_i     = Bi.*( mag_pi(:,ones1).*cos_k-mag_pk(:,ones1).*cos_i./delta);
        C_i     = mag_pi(:,ones1).^2-mag_pk(:,ones1).^2 ./delta;

        A=[A_k;A_k]; B=2*[B_i;B_k]; C=[C_i;C_k];

        % calling routine for determining solutions to quadratic equation
        [g1,g2] = quadrtic(A,B,C); size_g1 = size(g1);
        cbdb=[[[g1';g2'],bmag(:,j)];ones(2,size_g1(1)+1)];
        [abvblw, state(j)] = sectbnds(cbdb,0);
        bmag(:,j) = abvblw(1:2*m);
    end

    z=find(bmag(:,j)~=myeps & bmag(:,j)~=1/myeps & ...
        bmag(:,j)~=-248 & bmag(:,j)~=248 & ...
        bmag(:,j)~=-302 & bmag(:,j)~=302);
    if length(z), bmag(z,j)=bmag(z,j)*uP(nom(1),pct)*uC(nom(2),cct); end
    bmag(m+1:2*m,j)=(ones1')*(1/myeps);

    frac=j/length(pbds);
    set(info(2),'xdata',[0,frac,frac,0]);
    set(info(3),'string',[int2str(floor(100*frac)),'%']);
    drawnow;

    j=j+1;
end
close(info(1));

bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302) = ...
    20*log10(bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302));
bnd=[bmag; w(pbds); 7*ones(1,length(pbds))];

mesgbnds(w,pbds,state,state2,ph_r,bnd);

[jk,t]=sort(w(pbds));
if nargout==0
    plotbnds(bnd(:,t),[],180/pi*ph_r);
    title('SISOBND7 Bounds'),xlabel('X: Phase (degrees)  Y: Magnitude (dB)')
else
    bdb=bnd(:,t);
end
