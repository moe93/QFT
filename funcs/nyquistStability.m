function [ struct_out ] = nyquistStability( T_OL, DEBUG )
%
% Function to calculate closed-loop system stability according to
%   method introduced by Mario Garcia-Sanz (2016)
% 
% - Code free-download at:
% http://cesc.case.edu/Stability_Nyquist_GarciaSanz.htm
% And also included in the QFTCT, Controller design window
% -----------------------------------------------------------------
%
%   INPUTS:-
%       - T_OL      : Open-loop transfer function, i.e. L(s) = P(s) C(s)
%       - DEBUG     : If true, prints extra variables to help debug
%
%   RETURN:-
%       Returned variables are stored within the struct (u)
%       - zc        : = N + num_p_RHP
%       - N         : Total number of encirclements, N = Na + Nb + Nc + Nd
%       - num_p_RHP : Number of poles on the RHP in the s-plane
%       - Na        : Number of encirclements uncer case (Na)
%       - Nb        : Number of encirclements uncer case (Nb)
%       - Nc        : Number of encirclements uncer case (Nc)
%       - Nd        : Number of encirclements uncer case (Nd)
%       - zpCancel  : RHP zero-pole cancelations. "0"=no, "1"=yes
%       - k         : = [0, ±1, ±2, ...]depending on phase at w0 and w1
%       - sigma     : Describes sign of gain of T_OL ( = 0 if gain >= 0, 1
%                       otherwise)
%       - alpha     : = [0, +1, +2, ...]depending on phase at w0, w1, and
%                       w_inf
%       - gamma     : = [-2, 0, +2]depending on type and number of poles
%                       and zeros
%

arguments
    T_OL                    tf
    DEBUG                   {mustBeNumericOrLogical} = false
end

L = T_OL;                   % Set T_OL to L
% 1. Initial values
% -----------------
tolPh = 1e-3;               % tolerance phase
tolPh1 = 0.002;             % tolerance phase
tolPh2 = 5;                 % tolerance phase
tolM = 0.05;                % tolerance magnitude
tolA = 0.05;                % tolerance for change around 0.
wmin = 1e-16;               % lowest freq.
zpCancel = 0;               % RHP zero-pole cancelations. "0"=no, "1"=yes
z_LHP = []; z_RHP = []; z_0 = []; z_i = [];
p_LHP = []; p_RHP = []; p_0 = []; p_i = [];

% 2. gain, zeroes, poles, delay
% -----------------------------
dc_gain = abs(dcgain(L));
[~, ~, kkk]=tf2zpk(L.num{1},L.den{1});
zzz_roots = roots([L.num{1}]);
ppp_roots = roots([L.den{1}]);
orderNum = length(zzz_roots);
orderDen = length(ppp_roots);

timeDelay = L.iodelay;
zzzRe = real(zzz_roots);
pppRe = real(ppp_roots);
zzzMag = abs(zzz_roots);
pppMag = abs(ppp_roots);

[row_z_LHP,column_z_LHP] = find(zzzRe<0);
num_z_LHP = length(row_z_LHP);
for ii=1:num_z_LHP
    z_LHP(ii,1) = zzz_roots(row_z_LHP(ii),column_z_LHP(ii));
end

[row_z_RHP,column_z_RHP] = find(zzzRe>0);
num_z_RHP = length(row_z_RHP);
for ii=1:num_z_RHP
    z_RHP(ii,1) = zzz_roots(row_z_RHP(ii),column_z_RHP(ii));
end

[row_z_0,column_z_0] = find(zzzMag==0);
num_z_0 = length(row_z_0);
for ii=1:num_z_0
    z_0(ii,1) = zzz_roots(row_z_0(ii),column_z_0(ii));
end

[row_z_i,column_z_i] = find(zzzRe==0 & zzzMag~=0);
num_z_i = length(row_z_i);
for ii=1:num_z_i
    z_i(ii,1) = zzz_roots(row_z_i(ii),column_z_i(ii));
end
[row_p_LHP,column_p_LHP] = find(pppRe<0);
num_p_LHP = length(row_p_LHP);
for ii=1:num_p_LHP
    p_LHP(ii,1) = ppp_roots(row_p_LHP(ii),column_p_LHP(ii));
end

[row_p_RHP,column_p_RHP] = find(pppRe>0);
num_p_RHP = length(row_p_RHP);
for ii=1:num_p_RHP
    p_RHP(ii,1) = ppp_roots(row_p_RHP(ii),column_p_RHP(ii));
end

[row_p_0,column_p_0] = find(pppMag==0);
num_p_0 = length(row_p_0);
for ii=1:num_p_0
    p_0(ii,1) = ppp_roots(row_p_0(ii),column_p_0(ii));
end

[row_p_i,column_p_i] = find(pppRe==0 & pppMag~=0);
num_p_i = length(row_p_i);
for ii=1:num_p_i
    p_i(ii,1) = ppp_roots(row_p_i(ii),column_p_i(ii));
end

% 3. RHP zero-pole cancelation
% ---------------------------- (Sec.3-4. Rule 1)
cc = 0;
dd = [];
if( ~isempty(z_RHP) && ~isempty(p_RHP) )
    for pp=1:num_p_RHP
        dd = find(p_RHP(pp)==z_RHP);
        if ~isempty(dd)
            cc=1;
            break;
        end
    end
end
if (num_z_0>0 && num_p_0>0) || cc==1
    zpCancel = 1;
end

% 4. mag, pha, ww
% ---------------
% LN
[LNnum,LNden] = zp2tf([z_LHP;z_RHP],[p_LHP;p_RHP],kkk);
LN = tf(LNnum,LNden);
% [magLN0,phaseLN0] = bode(LN,0);
% L
[~, phaseL0] = bode( L, 1e-16 ); % phase at w=0
phaseL0 = round(phaseL0*100)/100; % Protection numerical accuracy
[~, pha2, ww2] = nichols(L);
ww1 = logspace(log10(ww2(1)),log10(ww2(end)),5000);
[mag1, pha1] = nichols( L, ww1 );
nn = length(pha1);
mag = [];
pha = [];
ww = [];
for jj=1:nn
    mag( jj )   = mag1( 1, 1, jj );
    pha( jj )   = pha1( 1, 1, jj );
     ww( jj )   =  ww1(    jj    );
end
% indPh2 = find( ~isequal(phaseL0, pha2), 1, 'first' ); % the first that is "~="
indPh2          = find( abs(phaseL0 - pha2) > tolPh, 1, 'first' ); % the first that is "~="
phaseL1         = pha2(indPh2); % phase at w=0+
leftRightAt0    = sign(phaseL0-phaseL1);

% 5. Find crosses at -900, -540, -180, +180 etc and mag>1
% -------------------------------------------------------
pha_neg360_0 = pha;
for jj=1:nn
    if pha_neg360_0(jj)>tolPh
        d1 = ceil(pha_neg360_0(jj)/360);
        pha_neg360_0(jj) = pha_neg360_0(jj) - 360*d1;
        
    elseif pha_neg360_0(jj)<=-360-tolPh
        d1 = floor(-pha_neg360_0(jj)/360);
        pha_neg360_0(jj) = pha_neg360_0(jj) + 360*d1;
    end
end
changeAround0 = pha_neg360_0 + 180;
indSignChange = [];
for jj=1:nn-1
    if (sign(changeAround0(jj+1))~=sign(changeAround0(jj))) & abs(pha_neg360_0(jj+1)-pha_neg360_0(jj))<(180-tolPh1) & abs(pha_neg360_0(jj+1)-pha_neg360_0(jj))>tolPh1
        indSignChange = [indSignChange jj];
    end
end

mm = length(indSignChange);
ind_kk_180 = [];
for jj=1:mm
    if( mag(indSignChange(jj)+1) > 1 && mag(indSignChange(jj)) > 1 )
        ind_kk_180 = [ind_kk_180 indSignChange(jj)];
    end
end

if ~isempty(ind_kk_180)
    n_ind_kk_180 = length(ind_kk_180);
end

% 6. Na. Fig.3.7(a)
% -----------------
Na = 0;
Na_1 = 0;
nMaxPha = length(pha);
if ~isempty(ind_kk_180)
    for jj=1:n_ind_kk_180
        Na_1(jj) = 0;
        if( mag(ind_kk_180(jj)) > (1+tolM) ) % if greater than 0 dB
            if (ind_kk_180(jj)+1)<nMaxPha % Protection
                pha21 = pha(ind_kk_180(jj)+1)-pha(ind_kk_180(jj));
                if abs(pha21)<tolPh
                    pha21 = 0;
                end
                ss = 1;
                while pha21==0
                    ss = ss+1;
                    if (nMaxPha-ind_kk_180(jj))>ss
                        pha21 = pha(ind_kk_180(jj)+ss)...
                            -pha(ind_kk_180(jj));
                        if abs(pha21)<tolPh
                            pha21 = 0;
                        end
                    else
                        pha21 = 0;

                        return;
                    end
                end
                if pha21>0
                    Na_1(jj) = -2; % Fig.3.7(a) to the right
                elseif pha21<0
                    Na_1(jj) = +2; % Fig.3.7(a) to the left
                else
                    Na_1(jj) = 0; % Fig.3.7(a) in axis
                end
            end
        end
    end
end
Na = sum(Na_1);
% if Na > 0
%     Na = ceil(Na/2); % Adjust Na to be an integer multiple of 2
% end

% 7. Nb. Fig.3.7(b)
% -----------------
Nb = 0;
if isfinite(dc_gain) % dc_gain is finite
    k_at_w0 = ((phaseL0/(-180))-1)/2; % at -900,-540,-180,+180,etc
    if abs(k_at_w0-round(k_at_w0))<tolPh
        if dc_gain>1 % dc_gain is >1
            pha21 = pha(2)-pha(1);
            ss = 2;
            while pha21==0
                ss = ss+1;
                if nMaxPha>ss
                    pha21 = pha(ss)-pha(1);
                else
                    pha21 = 0;
                    return;
                end
            end
            if pha21>0
                Nb = -1; % Fig.3.7(b) to the right
            elseif pha21<0
                Nb = +1; % Fig.3.7(b) to the left
            else
                Nb = 0; % Fig.3.7(b) in axis
            end
        else % dc_gain is <1
            Nb = 0; % Fig.3.7(b)
        end
    else
        Nb = 0; % Fig.3.7(b)
    end
end

% 8. Nc. Fig.3.7(c)
% -----------------
Nc = 0;

if( mag(end) > 1 )                      % mag( w=inf ) > 1
    k_at_wInf = ((pha(end)/(-180))-1)/2;
    if( abs(k_at_wInf-floor(k_at_wInf)) < 0.001 )
        % at -900,-540,-180,+180,etc
        pha21 = pha(end)-pha(end-1);
        ss = 1;
        while( pha21 == 0 )
            ss = ss+1;
            if( nMaxPha > ss )
                pha21 = pha(end)-pha(end-ss);
            else
                pha21 = 0;
                return;
            end
        end
        if( pha21 > 0 )
            Nc = -1;        % Fig.3.7(c) to the right
        elseif( pha21 < 0 )
            Nc = +1;        % Fig.3.7(c) to the left
        else
            Nc = 0;         % Fig.3.7(c) in axis
        end
    else
        Nc = 0; % Fig.3.7(c)
    end
end

% 9. Nd. Fig.3.7(d)
% -----------------
Nd = 0; k = 0; sigma = 0; alpha = 0; gamma = 0; % Initialization
if num_p_0 > 0
    if isfinite(dc_gain) % dc_gain is finite
        Nd = 0; % Fig.3.7(d)
    else % dc_gain is infinite
        % -- k --
        if( phaseL1 > -180 & phaseL1 < 90 )         % Case [a]
            k = -1;                                 %   Eq.(3.8)
        elseif( phaseL1 > 90 & phaseL1 < 180 )
            if( leftRightAt0 == 1 )                 % to the left. Case [c]
                k = 0;                              %   Eq.(3.10)
            elseif( leftRightAt0 == -1 )            % to the right. Case [b]
                k = -1;                             %   Eq.(3.9)
            end
        elseif( phaseL1 <= -180 )                   % to the left and right. Case [d]
            k = -ceil((phaseL1+180)/360);           %   Eq.(3.11)
        elseif( phaseL1 >= 180 )                    % to the left and right. Case [e]
            k = -ceil((phaseL1-180)/360);           %   Eq.(3.12)
        end
        % -- sigma --
        dcgainLN = dcgain(LN);

        if( dcgainLN >= 0 )                         % Case [a]
            sigma = 0;                               %   Eq.(3.13)
        else                                        % Case [b]
            sigma = 1;                               %   Eq.(3.14)
        end
        % -- gamma --
        numZP = num_z_RHP + num_p_LHP - num_z_LHP - num_p_RHP;
        gamma = 2 * numZP/max(abs(numZP),1); % Eq.(3.20)
        
        % -- alpha --
        if( phaseL0 <= 0 )
            if( leftRightAt0 == 1 )                 % to the left. Case [a]
                alpha = -ceil((phaseL0+90*(num_p_0))/360);
                                                    %   Eq.(3.15)
            elseif( leftRightAt0 == -1 )            % to the right. Case [b]
                alpha = -ceil((phaseL0+90*(num_p_0-2))/360);
                                                    %   Eq.(3.16)
            end
        elseif( phaseL0 > 0 )
            if( leftRightAt0 == 1 )                 % to the left. Case [c]
                alpha = floor((phaseL0+90*(num_p_0-2))/360);
                                                    %   Eq.(3.17)
            elseif( leftRightAt0 == -1 )            % to the right. Case [d]
                alpha = -floor((phaseL0-90*(num_p_0))/360);
                                                    %   Eq.(3.18)
            end
        end
        if( (  orderNum == orderDen  ) && ...
            ( (pha(end) < 180+tolPh2 && pha(end) >180-tolPh2) | ...
              (pha(end) <     tolPh2 && pha(end) >   -tolPh2)) ) % Case [e]
            alpha = 0; % Eq.(3.19)
        end
        Nd = 2*(k+1) + sigma + alpha*gamma; % Eq.(3.6)
    end
end

% 10. N sum
% ---------
N = Na + Nb + Nc + Nd; % Eq.(3.5)

% 11. Zc sum
% ---------- (Sec.3-4. Rule 2)
zc = N + num_p_RHP;

% The closed-loop system is stable if:
% - Rule 1: “zpCancel=0”. RHP zero-pole cancelations "0"=no,"1"=yes
% - Rule 2: “zc = 0”. Being zc = N + num_p_RHP
% ---------------------------------------------



% --- If DEBUG is flagged, print the following
if( DEBUG )
    disp( '======================================' );
    disp( '=========== DEBUG __ START ===========' );
    disp( '======================================' );
    if( exist('leftRightAt0', 'var') )
        fprintf( "\tDelta      = %6.3f\n",  leftRightAt0        );  % +ve if Lw goes initially to the left
    else
        fprintf( "\tDelta      = D.N.E\n"                       );
    end
    if( exist('dcgainLN', 'var') )
        fprintf( "\tDC Gain    = %6.3f\n",  dcgainLN            );
    else
        fprintf( "\tDC Gain    = D.N.E\n"                       );
    end
    if( exist('phaseL0', 'var') )
        fprintf( "\t∠ L(w0)    = %6.3f\n",  phaseL0             );
    else
        fprintf( "\t∠ L(w0)    = D.N.E\n"                       );
    end
    if( exist('phaseL1', 'var') )
        fprintf( "\t∠ L(w1)    = %6.3f\n",  phaseL1             );
    else
        fprintf( "\t∠ L(w1)    = D.N.E\n"                       );
    end
    if( exist('numZP', 'var') )
        fprintf( "\tnumZP      = %6.3f\n",  numZP               );
        fprintf( "\tmax(numZP) = %6.3f\n",  max(abs(numZP), 1)  );
    else
        fprintf( "\tnumZP      = D.N.E\n"                       );
        fprintf( "\tmax(numZP) = D.N.E\n"                       );
    end
    disp( '======================================' );
    disp( '============ DEBUG __ END ============' );
    disp( '======================================' );
end

% --- Aggregate output
value       = [ zc, N, num_p_RHP, Na, Nb, Nc, Nd, ...
                zpCancel, k, sigma, alpha, gamma].';
variable    = [ "zc", "N"   , "p_RHP", "Na" , "Nb"  , ...
                "Nc", "Nd"  , "zpCancel"    , "k"   , ...
                "sigma"     , "alpha"       , "gamma" ].';

% --- If no output is requested, print the values we computed
if( ~nargout )
    disp( table( variable, value ) )

% --- Otherwise, generate struct and output it
else
    namesoftags = cellstr(variable);
    for i = 1:length( namesoftags )
        struct_out.(namesoftags{i}) = value(i);
    end
end

end
