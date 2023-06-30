% QFTEX1 Main Example.

% Author: Y. Chait
% 10/10/94
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

%% Setup environment
clear all; close all; clc;
set( groot, 'defaultLineLineWidth', 1.5 );	% Set default line width of plots
set( 0, 'DefaultAxesFontSize', 12 );        % Set a default axes font size
% Change default interpreter (affects text, labels, etc...)
set( 0, 'DefaultTextInterpreter', 'latex' );
set( 0, 'DefaultLegendInterpreter', 'latex' );
set( 0, 'DefaultAxesTickLabelInterpreter', 'latex' );
format compact;
fontsize = 12;

% Flags
CNTR = 1;                                   % Figure handle counter
PLOT = true;                                %#ok<NASGU> If true, plot figures!
PLOT = false;                               % COMMENT OUT TO PRINT FIGURES
PRNT = true;                                %#ok<NASGU>
% PRNT = false;                               % COMMENT OUT TO PRINT FIGURES

%% Add folders/files to path
% Get current path to working directory and split
pathParts = strsplit( pwd, filesep );
% Go up one level and generate new path
src = fullfile( pathParts{1:end-1} );

% If on a UNIX machine (i.e. macOS, Ubuntu, etc...), fix path since
% strsplit() removes the leading '/' from the path.
if( isunix )
    src = [ filesep src ];
end

% Add QFT2 to path
addpath( genpath(src) );

%% [INFO] Strings
ACK = 'COMPLETED\n\n';

%%
% Example 2.1 from Mario Garcia-Sanz's:
%   "Robust Control Engineering: Practical QFT Solutions"
%
%
%	Consider a continuous-time, siso, negative unity feedback system
%
%                                        |V(s)                |D(s)
%          ----  +           ----        v        ----        v
%   ----> |F(s)|--> O ----> |G(s)| ----> O ----> |P(s)| ----> O ------>
%   R(s)   ----     ^ -      ----                 ----        |  Y(s)
%                   |                  ----                   |
%                   ----------------- |H(s)| -----------------
%                                      ----
%
%
%       The plant P(s) has parametric a uncertainty model:
%	                                 KT                    
%	   P(s)  =  ------------------------------------------- 
%	            (J.La)s^2 + (J.Ra + b.La)s + (b.Ra + KT.Ke)     
%
%   Where
%       J  in [ J_min,  J_max], b  in [ b_min,  b_max]
%       La in [La_min, La_max], Ra in [Ra_min, Ra_max]
%       KT in [KT_min, KT_max],         Ke = KT
%
%
%	The performance specifications are: design a controller
%       G(s) such that it achieves
%
%  1) Type 1: Stability specification
%       ==> Corresponds to sisobnds( 1, ... )
%
%                   |  P(jw)G(jw)  |
%       |T_1(jw)| = |--------------| <= del_1(w) == W_s == 1.46
%                   |1 + P(jw)G(jw)|
%
%       w in [1, 10000] rad/s
%
%  2) Type 3: Sensitivity OR Disturbance at plant output specification 
%       ==> Corresponds to sisobnds( 2, ... )
%
%                   |       1      |                  (s/a_d)
%       |T_3(jw)| = |--------------| <= del_3(w) == -----------
%                   |1 + P(jw)G(jw)|                (s/a_d) + 1
%
%       w in [1, 1000] rad/s ; a_d = 1,000
%
%  3) Type 6: Reference tracking specification
%       ==> Corresponds to sisobnds( 7, ... )
%
%                                 |        P(jw)G(jw)  |
%       del_6_lo(w) < |T_3(jw)| = |F(jw) --------------| <= del_6_hi(w)
%                                 |      1 + P(jw)G(jw)|
%
%       w in [1, 5000] rad/s
%
%   Where
%
%                     |        1        |
%       del_6_lo(w) = |-----------------| ; a_L = 1000
%                     | ((s/a_L) + 1)^2 |
%
%                     |       ((s/a_U) + 1)       |
%       del_6_hi(w) = |---------------------------| ; a_U = 1000
%                     | (s/w_n)^2 + (2zs/w_n) + 1 |
%
%       z = 0.8 , w_n = 1.25a_u/z
%
%

%% Step 1: Plant Modeling & Uncertainty

% --- Plant parameters
%   min_    : Minimum value
%   max_    : Maximum value
%   grid_   : Gridding
%
min_k    = 0.2; max_k    = 2.0 ; grid_k    = 2;
min_z    = 0.5; max_z    = 0.75; grid_z    = 2;
min_p    = 1.0; max_p    = 10.0; grid_p    = 2;
min_wn   = 5.0; max_wn   = 6.0 ; grid_wn   = 2;
min_zeta = 0.8; max_zeta = 0.9 ; grid_zeta = 2;


% --- Gridding
%   ***NOTE: Can grid using logspace() or linspace()
%   _g  : Gridded variable
%
k_g    = logspace( log10(min_k)    , log10(max_k)      , grid_k    );
z_g    = logspace( log10(min_z)    , log10(max_z)      , grid_z    );
p_g    = logspace( log10(min_p)    , log10(max_p)      , grid_p    );
wn_g   = logspace( log10(min_wn)   , log10(max_wn)     , grid_wn   );
zeta_g = logspace( log10(min_zeta) , log10(max_zeta)   , grid_zeta );


% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 900 ) == SISO with 900 TFs
%
n_Plants = grid_k*grid_z*grid_p*grid_wn*grid_zeta;   % Number of plants
P = tf( zeros(1,1,n_Plants) );                      % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

syms s;
NDX = 1;                                            % Plant counter
for var1 = 1:grid_k                                 % Loop over k
    k = k_g( var1 );                                % ....

    for var2 = 1:grid_z                             % Loop over z
        z = z_g( var2 );                            % ....

        for var3 = 1:grid_p                         % Loop over p
            p = p_g( var3 );                        % ...

            for var4 = 1:grid_wn                    % Loop over wn
                wn = wn_g( var4 );                  % ...

                for var5 = 1:grid_zeta              % Loop over zeta
                    zeta = zeta_g( var5 );          % ...

                    % --- Here we create the plant TF
                    num = sym2poly( k*(s/z + 1) );
                    den = sym2poly( s*(s/p + 1)*((s/wn)^2 + (2*zeta/wn)*s + 1) );
                    P(1, 1, NDX) = tf( num, den );  % Transfer Function
                    NDX = NDX + 1;                  % Incerement counter
                end
            end
        end
    end
end
clear s;

% [INFO] ...
fprintf( ACK );

%% Step 2: The Nominal Plant
k_0     = 0.2;
z_0     = 0.5;
p_0     = 1.0;
wn_0    = 5.0;
zeta_0  = 0.8;

% [INFO] ...
fprintf( 'Step 2:' );
fprintf( '\tComputing nominal plant...' );

% --- Generate nominal plant TF
%   Any one of the models above can be used as the nominal plant.
%   We just happened to chose this one.
%
syms s;
num = sym2poly( k_0*(s/z_0 + 1) );
den = sym2poly( s*(s/p_0 + 1)*((s/wn_0)^2 + (2*zeta_0/wn_0)*s + 1) );
clear s;
P_0(1, 1, 1) = tf( num, den );          % Nominal Transfer Function

% --- Alternatively, just use the following
P_0(1, 1, 1) = P( :, :, 1 );            % Nominal Transfer Function
% --- Define nominal plant case
nompt = 1;

% % --- Append to the end of the gridded plants
% P( 1, 1, end+1 ) = P_0;
% % --- Define nominal plant case
% nompt = length( P );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = logspace( log10(0.001), log10(300), 1024);
bode( P_0, w ); grid on;

if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( P_0, w ); grid on;
    make_nice_plot();
end

% --- Plot root locus
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    rlocus( P_0 );
    title('Root Locus of Plant (under Proportional Control)')
    make_nice_plot();
end

%% Step 3: QFT Template

% [INFO] ...
fprintf( 'Step 3:' );
fprintf( '\tPlotting QFT templates...' );

% --- Working frequencies
w = [ 0.001 0.01 0.1 0.5 1 4 10 30 150 300 ];

if( PLOT )
    % --- Plot QFT templates
    plottmpl( w, P, nompt );
    
    % --- Change legend position
    hLegend = findobj( gcf, 'Type', 'Legend' ); % Get legend property
    set( hLegend, 'location', 'southeast' );    % Access and change location
    
    % --- Change plot limits
%     xmin = -25; xmax = 10; dx = 5;
%     xlim( [xmin xmax] );
%     xticks( xmin:dx:xmax )
    title( 'Plant Templates' )
    
    % --- Beautify plot
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

%% Step 4: Define Stability Specifications
% --- The stability margins, gain and phase, are specified here
%

% L(s) = G(s)P(s) where G(s) == control, P(s) == plant
%
%                L(s)        G(s)P(s)
% Then H(s) = --------- = --------------
%              1 + L(s)    1 + G(s)P(s)  
%

% [INFO] ...
fprintf( 'Step 4:' );
fprintf( '\tDefining stability specifications...' );

% --- Type 1
% Frequencies of interest
omega_1 = [ 0.001 0.01 0.1 0.5 1 4 10 30 150 300 ];
% Restriction
W_s         = 1.2;
del_1       = W_s;
PM          = 180 - 2*(180/pi)*acos(0.5/W_s);        % In deg
GM          = 20*log10( 1+1/W_s );                   % In dB

% [INFO] ...
fprintf( ACK );

%% Step 5: Define Performance Specifications

% [INFO] ...
fprintf( 'Step 5:' );
fprintf( '\tDefining performance specifications...' );

% --- Type 3: Sensitivity or ouptut disturbance rejection specification
% Frequencies of interest
omega_3 = [ 0.001 0.01 0.1 0.5 ];

% Restriction
num     = [ 0.5 2 2.8 0 ];
den     = [ 0.5 2 2.8 1 ];
del_3   = tf( num, den );


% --- Type 6: Reference tracking specification
% Frequencies of interest
omega_6 = [ 0.001 0.01 0.1 0.5 1 4 10 ];

% Restriction
% Upper bound
num = [ 2.8571 1 ];
den = [ 0.6667 2.3333 1 ];
del_6_hi = tf( num, den );
% Lower bound
num = 1;
den = [ 0.5 2 2.8 1 ];
del_6_lo = tf( num, den );
% Tracking weight
del_6 = [ del_6_hi  ;
          del_6_lo ];

% [INFO] ...
fprintf( ACK );

%% Step 6: Calculate Staibility QFT Bounds

% --- Example 2.1 continued (Pg. 36)
%   - Type 1: Stability specification
%       > Corresponds to sisobnds( 1, ... )
%   - Type 3    //
%       > Corresponds to sisobnds( 2, ... )
%   - Type 6    //
%       > Corresponds to sisobnds( 7, ... )
%

% --------------------------------------------------
% ----      Type 1: Stability specification     ----
% --------------------------------------------------
spec = 1;

% [INFO] ...
fprintf( 'Step 6:' );
fprintf( '\tCalculating stability QFT bounds\n' );
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb1 = sisobnds( spec, omega_1, del_1, P, [], nompt );
% R = 0; bdb1 = sisobnds( spec, omega_1, del_1, P, R, nompt );

if( PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );
    
    % --- Plot bounds
    plotbnds( bdb1 );
    title( 'Robust Stability Bounds' );
%     xlim( [-360 0] ); ylim( [-10 30] );
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

%% Step 7: Calculate Performance QFT Bounds

% -------------------------------------------
% ---- Type 3: Sensitivity specification ----
% -------------------------------------------
spec = 2;

% [INFO] ...
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb2 = sisobnds( spec, omega_3, del_3, P, [], nompt );

if( PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );
    
    % --- Plot bounds
    plotbnds(bdb7);
    title( 'Sensitivity Reduction Bounds' );
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

% --------------------------------------------------
% ---- Type 6: Reference tracking specification ----
% --------------------------------------------------
spec = 7;

% [INFO] ...
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb7 = sisobnds( spec, omega_6, del_6, P );

if( PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );
    
    % --- Plot bounds
    plotbnds(bdb7);
    title( 'Robust Tracking Bounds' );
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );


%% Step 8: Intersection of QFT Bounds and Compatibility

% [INFO] ...
fprintf( 'Step 8:' );
fprintf( '\tGrouping bounds...' );

% --- Grouping bounds
bdb = grpbnds( bdb1, bdb2, bdb7 );
if( PLOT )
    plotbnds( bdb );
    title('All Bounds');
end

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

% --- Find bound intersections
ubdb = sectbnds( bdb );
if( PLOT )
    plotbnds( ubdb );
    title('Intersection of Bounds');
end

% [INFO] ...
fprintf( ACK );

%% Step 9: Synthesize Feedback Controller G(s)

% [INFO] ...
fprintf( 'Step 9:' );
fprintf( '\tSynthesize G(s)...' );

% --- Directory where QFT generated controllers are stored
src = './controllerDesigns/';
% --- Controller file, G(s)
G_file  = [ src 'airplaneFlightController.shp' ];
if( isfile(G_file) )
    G = getqft( G_file );
else
    syms s;
    num = 5 .* sym2poly( (s/4 + 1)*(s/13 + 1) );    % Get coefficients
    den = 1 .* sym2poly( (s/350 + 1)^2 );           % ...
    clear s;
    
    % Construct controller TF
    G = tf( num, den );
end

% Define a frequency array for loop shaping
wl = logspace( log10(0.001), log10(300), 2048 );
L0 = P( 1, 1, nompt );
L0.ioDelay = 0; % no delay
lpshape( wl, ubdb, L0, G );

% [INFO] ...
fprintf( ACK );


%% Step 10: Synthesize Prefitler F(s)

% [INFO] ...
fprintf( 'Step 10:' );
fprintf( '\tSynthesize F(s)...' );

% --- Directory where QFT generated controllers are stored
src = './controllerDesigns/';
% --- Pre-filter file, F(s)
F_file  = [ src 'airplaneFlightController.fsh' ];
if( isfile(F_file) )
    F = getqft( F_file );
else
    syms s;
    num = 1;                                        % Get coefficients
    den = 1 .* sym2poly( (s/3.5 + 1) );             % ...
    clear s;
    
    % Construct controller TF
    F = tf( num, den );
end

pfshape( 7, wl, del_6, P, [], G, [], F );

% [INFO] ...
fprintf( ACK );

%% Step 11-13: ANALYSIS

disp(' ')
disp('chksiso(1,wl,del_1,P,R,G); %margins spec')
% chksiso( 1, wl, del_1, P, [], G );
chksiso( 1, wl, del_1, P, [], G, [], F );

disp(' ')
disp('chksiso(2,wl,del_3,P,R,G); %Sensitivity reduction spec')
ind = find(wl <= 0.5);
% chksiso( 2, wl(ind), del_3, P, [], G );
chksiso( 2, wl(ind), del_3, P, [], G, [], F );

disp(' ')
disp('chksiso(7,wl,W3,P,R,G); %input disturbance rejection spec')
drawnow
% chksiso(7,wl(ind),W7,P,[],G);
chksiso( 7, wl, del_6, P, [], G, [], F );


%% Feedback and feedforward
% Recall, G_f(s) = -P_c(s)^-1*M(s)*V(s)
%   Where,
%
%       1) G_f(s): Feedforward element
%       2) P_c(s): A central plant within the uncertainty
%       3)   M(s): Dynamics of the disturbance over the plant output
%       3)   V(s): A LP filter with high frequency poles
%

% --- Define P_c(s)
k_c     = 1.10 ;
z_c     = 0.625;
p_c     = 5.50 ;
wn_c    = 5.50 ;
zeta_c  = 0.85 ;

syms s;
num = sym2poly( k_c*(s/z_c + 1) );
den = sym2poly( s*(s/p_c + 1)*((s/wn_c)^2 + (2*zeta_c/wn_c)*s + 1) );
clear s;
P_c(1, 1, 1) = tf( num, den );          % Central plant within uncertainty

% --- Define M(s)
M = tf( 0.2, conv([1 0], [1 1]) );

% --- Define V(s)
V = tf( 1, [1/5.5 1] );

% --- Construct G_f(s)
G_f = minreal( -P_c^-1*M*V );

% --- Update G(s) to consider contribution of G_f(s)
syms s;
num = 60 .* sym2poly( (s/6 + 1)*(s/10 + 1)*(s/150 + 1) );
den = sym2poly( (s/70+1)*(s/200+1)^2 );
clear s;
G_updated = tf( num, den );         % Updated controller

%% Test the feedforward bound generator
% Recall, the disturbance rejection with feedforward element is defined as
%
%   | y(s) |   | M(s) + P(s)G_f(s) |
%   | ---- | = | ----------------- | <= δ_d(𝜔) , 𝜔 𝜖 Ω_d
%   | d(s) |   |    1 + P(s)G(s)   |
%
%   bdb = genbnds( ptype, w, Ws, A, B, C, D, Pnom, phs );
%
%   | y(s) |   |  A(s) + B(s)G(s)  |
%   | ---- | = | ----------------- | <= Ws
%   | d(s) |   |  C(s) + D(s)G(s)  |
%

% margins (1,1): g2=0
disp(' ')
disp( 'bdb12=genbnds(12,w,W11,a,b,c,d,P(1,1,1));' );

% --- Working frequencies
w = [ 0.001 0.01 0.1 0.5 1 4 10 ];

% --- Desired disturbance rejection specification
W11 = 0.014;

% --- Matrices
a = M;
b = P*G_f;
c = 1;
d = P;

% --- Bound generation (ptype = 12: Sensitivity or output disturbance
%                                   rejection with feedforward)
bdb12 = genbnds( 12, w, W11, a, b, c, d, P_0 );

% % --- Plot bounds
plotbnds(bdb12);
title('Disturbance rejection specification with feedforward element Bounds');
make_nice_plot();

%% Step XYZ: Re-intersect QFT Bounds with new feedforward  bound

% [INFO] ...
fprintf( 'Step 8:' );
fprintf( '\tGrouping bounds...' );

% --- Grouping bounds
% bdb = grpbnds( bdb1, bdb12, bdb7 );
bdb = grpbnds( bdb1, bdb12);
if( ~PLOT )
    plotbnds( bdb );
    title('All Bounds');
end

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

% --- Find bound intersections
ubdb = sectbnds( bdb );
if( ~PLOT )
    plotbnds( ubdb );
    title('Intersection of Bounds');
end

% [INFO] ...
fprintf( ACK );

%% Misc

syms s;
num = -0.1818 .* sym2poly( (s/5.5)^2 + (2*0.85/5.5)*s + 1 );
den = sym2poly( (s/0.625 +1)*(s + 1) );
Gf = tf( num, den );
clear s;
