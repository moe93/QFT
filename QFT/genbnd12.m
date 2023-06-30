function bdb = genbnd12( w, wbd, W, M, PG_f, C, P, Pnom, ph_r, info, maxr )
% function bdb = genbnd12( w, wbd, W, A, B, C, D, Pnom, ph_r, info, maxr )
% Compute QFT bounds for the following closed-loop configuration
%
%                          |  A + BG_f  |
%                          | ---------- | <= WS
%                          |   C + DG   |
%
%     GENBNDS(10,W,WBD,WS,A,B,C,D,Pnom,PHS) computes bounds for the
%     above configuration.  W is the entire set of
%     frequencies and WBD is a subset of W designating bounds to compute.
%     WS is the performance specification. A, B, C and D are
%     the frequency response data sets, Pnom designates the nominal plant
%     for the configuration. PHS specifies at which phases (degrees)
%     to compute bounds.
%

% Author: Craig Borghesani
% 1/27/94
% Copyright (c) 2003, Terasoft, Inc.

% --- Re-assign variables so I don't have to re-write a bunch of code
A = M;
B = PG_f;
C = C;
D = P;

disp( 'Hey there, ptype = 12.' );
mLen    = length( ph_r );       % Number of phases to go through
pbds    = qsubset( wbd, w );
state2  = [];
myeps   = 1e-16;
% --- Unused variables
% final_state = ones( 1, length(wbd) );
% nbd     = length( w );          % Number of bounds

[ rma, cma ] = size( A ); [ rmc, cmc ] = size( C );
[ rmb, cmb ] = size( B ); [ rmd, cmd ] = size( D );
[ rmw,  ~  ] = size( W );
row = [ rma, rmb, rmc, rmd, rmw ];
maxr = max( row );

% preallocating matrix for more efficient programming
bmag = [ ones( mLen, length(pbds) )*myeps    ;
         ones( mLen, length(pbds) )*1/myeps ];

% --- Get phase (u_) and angle (v_)
uPnom   = abs( Pnom ); vPnom = qatan4( Pnom );
uA      = abs( A )   ; vA    = qatan4( A )   ;
uB      = abs( B )   ; vB    = qatan4( B )   ;
uC      = abs( C )   ; vC    = qatan4( C )   ;
uD      = abs( D )   ; vD    = qatan4( D )   ;

% declaring sizes of replicating matricies
if( repltest )
    u = ones( maxr, 1    );
    v = ones( 1   , mLen );
else
    u = ones( 1   , maxr );
    v = ones( 1   , mLen );
end

j = 1; act = 1; bct = 1; cct = 1; dct = 1;
while( j <= length(pbds) )

    phi = ph_r(u, :) - (vPnom(pbds(j)));
    
    % --- cma == Columns of Matrix A
    if( cma > 1 ); act = pbds(j); end
    if( cmb > 1 ); bct = pbds(j); end
    if( cmc > 1 ); cct = pbds(j); end
    if( cmd > 1 ); dct = pbds(j); end
    
    % --- Get magnitudes of matrices (m_) for frequency w_j
    mA = uA( :, act ); mB = uB( :, bct );
    mC = uC( :, cct ); mD = uD( :, dct );

    % --- Get phases     of matrices (p_) for frequency w_j
    pA = vA( :, act ); pB = vB( :, bct );
    pC = vC( :, cct ); pD = vD( :, dct );

    %%%%%% V5 code
    if( length(mA) == 1 ); mA = mA( u ); end
    if( length(mB) == 1 ); mB = mB( u ); end
    if( length(mC) == 1 ); mC = mC( u ); end
    if( length(mD) == 1 ); mD = mD( u ); end
    if( length(pA) == 1 ); pA = pA( u ); end
    if( length(pB) == 1 ); pB = pB( u ); end
    if( length(pC) == 1 ); pC = pC( u ); end
    if( length(pD) == 1 ); pD = pD( u ); end
    
    % --- Setup angles (phases)
    alpha = pB + phi;   %alpha = alpha( :, v ) - phi;
    beta  = pA + phi;   %beta  = beta( :, v ) - phi;
    psi   = pD + phi;   %psi   = psi( :, v ) - phi;
    % --------------------------------------------
    
    % --- Extract distrubance
    Ws = W( :, pbds(j) );
    
    % --- Setup quadratic equation ( i.e. a*g^2 + b*g + c >= 0 )
    % Note:
    %   |     M(jw)    |    = m     = mA
    %   | P(jw)G_f(jw) |    = pg_f  = mB
    %   |     C(jw)    |    = c     = mC
    %   |     P(jw)    |    = p     = mD
    %
    
    % --- Extract variables into easier to use names
    m    = mA( :, v );
    pg_f = mB( :, v );
    % c    = mC( :, v );
    p    = mD( :, v );

    % --- Place into quadratic coefficients
    a = p.^2;
    b = 2.*p.*cos(psi);
    c = -( m.^2 + pg_f.^2 + 2.*m.*pg_f.*cos(alpha-beta) ) ./ Ws(:,v).^2 + 1;
    % --------------------------------------------
    
    [ g1, g2 ] = quadrtic( a, b, c ); size_g1 = size(g1);
    cbdb = [ [[g1';g2'], bmag(:,j)] ;
             ones(2,size_g1(1)+1)  ];
    [ abvblw, state(j) ] = sectbnds( cbdb, 0 );
    bmag( :, j ) = abvblw( 1:2*mLen ) ;

    z = find( bmag(:,j) ~= myeps & bmag(:,j) ~= 1/myeps & ...
              bmag(:,j) ~= -248  & bmag(:,j) ~= 248     & ...
              bmag(:,j) ~= -302  & bmag(:,j) ~= 302 );
    if( ~isempty(z) )
        bmag(z,j) = bmag(z,j)*uPnom(pbds(j));
    end

    frac = j/length(pbds);
    set(info(2),'xdata',[0,frac,frac,0]);
    set(info(3),'string',[int2str(floor(100*frac)),'%']);
    drawnow;

    j=j+1;
end
close(info(1));

bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302) = ...
    20*log10(bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302));

bnd=[bmag; w(pbds); 10*ones(1,length(pbds))];

mesgbnds(w,pbds,state,state2,ph_r,bnd);

[jk,t]=sort(w(pbds));
if nargout==0
    plotbnds(bnd(:,t),[],180/pi*ph_r);
    title('GENBND12 Bounds'),xlabel('X: Phase (degrees)  Y: Magnitude (dB)')
else
    bdb=bnd(:,t);
end
