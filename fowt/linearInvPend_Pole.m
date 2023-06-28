% Linearized inverted pendulum QFT control
%   Pole - Full (does NOT assume M >> m )
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Jun.  5th, 2023
%
% CHANGELOG :
%   Jun.  5th, 2023
%       - Initial script
%

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

%% Read A, B, C, D matrices from linearized model
data_dir    = './data/';
name_mdl    = 'SS_linearizedInvertedPendulum.mat';
stateSpace  = load( [data_dir name_mdl ] );

% --- Get number of states
nStates = stateSpace.nx;
% --- Extract matrices from the one, big ABCD matrix
A_full = stateSpace.ABCD( 1:nStates      , 1:nStates     );
B_full = stateSpace.ABCD( 1:nStates      , nStates+1:end );
C_full = stateSpace.ABCD( nStates+1:end  , 1:nStates     );
D_full = stateSpace.ABCD( nStates+1:end  , nStates+1:end );

% --- In this case, we only care about the first 4 states
nStatesKeep = 4;
A = A_full( 1:nStatesKeep   , 1:nStatesKeep );
B = B_full( 1:nStatesKeep   , 1:end         );
C = C_full( 1:end           , 1:nStatesKeep );
D = D_full( 1:height(C)     , 1:end         );

% --- Generate state-space model
% States and inputs names
stateNames  = [ "x" "x_dot" "theta" "theta_dot" ];
inputNames  = [ "Fx" ];
outputNames = [ "x" "theta" ];
% State-space model
sys         = ss( A, B, C, D                , ...
                  'StateName' , stateNames  , ...
                  'InputName' , inputNames  , ...
                  'OutputName', outputNames );

% --- Generate TF from SS model
TF = tf( sys );

%% DOUBLE CHECK
% Generate SS model using analytical approach

M_0 = 2.0   ;                   % Mass of cart                  [  kg  ]
m_0 = 0.075 ;                   % Mass of rod                   [  kg  ]
g   = 9.8   ;                   % Gravitational acceleration    [m.s^-2]
h   = 0.5   ;                   % Rod length                    [  m   ]

% ------------------------------
A_theory = [ 0      1       0                   0   ;
             0      0   -m_0*g/M_0              0   ;
             0      0       0                   1   ;
             0      0   (M_0+m_0)*g/(M_0*h)     0 ] ;
% ------------------------------
B_theory = [ 0              ;
             1/M_0          ;
             0              ;
             -1/(M_0*h) ]   ;
% ------------------------------
C_theory = [ 1  0   0   0   ;
             0  0   1   0 ] ;
% ------------------------------
D_theory = [ 0  ;
             0 ];

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};

sys_theory = ss( A_theory, B_theory, C_theory, D_theory , ...
                 'StateName'    ,   stateNames          , ...
                 'InputName'    ,   inputNames          , ...
                 'OutputName'   ,   outputNames );

% --- Generate TF from SS model
TF_theory = tf( sys_theory );

%% Step 1: Plant Modeling & Uncertainty

% --- Plant parameters
%   min_    : Minimum value
%   max_    : Maximum value
%   grid_   : Gridding
%
min_M   = 1.50;     max_M   = 2.50;     grid_M  = 5;
min_m   = 0.05;     max_m   = 0.09;     grid_m  = 5;


% --- Gridding
%   ***NOTE: Can grid using logspace() or linspace()
%   _g  : Gridded variable
%
M_g = logspace( log10(min_M)    ,   log10(max_M)    ,   grid_M );
m_g = logspace( log10(min_m)    ,   log10(max_m)    ,   grid_m );


% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 300 ) == SISO with 300 TFs
%
n_Plants = grid_M*grid_m;                           % Number of plants
P_y1 = tf( zeros(1,1,n_Plants) );                   % Pre-allocate memory
P_y2 = tf( zeros(1,1,n_Plants) );                   % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

NDX = 1;                                            % Plant counter
for var1 = 1:grid_M                                 % Loop over M
    M = M_g( var1 );                                % ....

    for var2 = 1:grid_m                             % Loop over m
        m = m_g( var2 );                            % ....

        % --- Here we create the plant TF
        A_g = [ 0      1       0            0   ;
                0      0   -m*g/M           0   ;
                0      0       0            1   ;
                0      0   (M+m)*g/(M*h)    0 ] ;
        % ------------------------------
        B_g = [ 0       ;
                1/M     ;
                0       ;
               -1/(M*h) ];
        % ------------------------------
        C_g = [ 1 0 0 0  ;
                0 0 1 0 ];
        % ------------------------------
        D_g = [ 0  ;
                0 ];

        % --- Generate grided TF from grided SS model
        sys_g = ss( A_g, B_g, C_g, D_g              , ...
                    'StateName'    ,   stateNames   , ...
                    'InputName'    ,   inputNames   , ...
                    'OutputName'   ,   outputNames );
        TF_g = tf( sys_g );
        P_y1(:, :, NDX) = TF_g(1);      % Transfer Function
        P_y2(:, :, NDX) = TF_g(2);      % Transfer Function
        NDX = NDX + 1;                  % Incerement counter
    end
end

% --- Compare to pre-known form from MGS book
syms s;
P_11 = vpa( (-h*s^2+g) / (s^2*(-M_0*h*s^2+(M_0+m_0)*g)) );  % Cart/Input
P_21 = vpa( 1 / (-M_0*h*s^2+(M_0+m_0)*g) );                 % Pole/Input
clear s;

% [INFO] ...
fprintf( ACK );

%% Step 2: The Nominal Plant

% --- Re-state those variables just in case we forgot
% M_0 = 2000.0  ;                   % Mass of cart                  [  kg  ]
% m_0 = 375.00  ;                   % Msas of rod                   [  kg  ]
% b_0 = 0.2000  ;                   % Co-efficient of friction
% I_0 = 1126.95 ;                   % 2nd mass moment of inertia    [kg.m^2]
% g_0 = 9.80665 ;                   % Gravitational acceleration    [m.s^-2]
% l_0 = 6/2     ;                   % Half the rod length           [  m   ]

% [INFO] ...
fprintf( 'Step 2:' );
fprintf( '\tComputing nominal plant...' );

% --- Generate nominal plant TF
%   Any one of the models above can be used as the nominal plant.
%   We just happened to chose this one.
%
P_y10(1, 1, 1) = TF_theory(1);          % Nominal Transfer Function
P_y20(1, 1, 1) = TF_theory(2);          % Nominal Transfer Function

% --- Append to the end of the gridded plants
P_y1( 1, 1, end+1 ) = P_y10;
P_y2( 1, 1, end+1 ) = P_y20;

% --- Select which I/O pair we want to work with
P   = P_y2  ;                           % This corresponds to gridded x/Fx
P_0 = P_y20 ;                           % This corresponds to nominal x/Fx
SS_0 = ss(P_0);                         % Convert to SS for use in Simulink
% --- Define nominal plant case
nompt = length( P );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = logspace( -4, 3 );
[p0, theta0] = bode( P_0, w );

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
w =  [ 0.01 0.05 0.1 0.5 1 5 10 50 100 500 ];

if( PLOT )
    % --- Plot QFT templates
    plottmpl( w, P, nompt );
    
    % --- Change legend position
    hLegend = findobj( gcf, 'Type', 'Legend' ); % Get legend property
    set( hLegend, 'location', 'southeast' );    % Access and change location
    
    % --- Change plot limits
    xmin = -25; xmax = 10; dx = 5;
    xlim( [xmin xmax] );
    xticks( xmin:dx:xmax )
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
fprintf( '\tDefining stability specifications\n' );

% --- Type 1
% Frequencies of interest
omega_1 = [ 0.01 0.05 0.1 0.5 1 5 10 50 100 500 ];
% Restriction
W_s         = 1.08;
del_1       = W_s;
PM          = 180 -2*(180/pi)*acos(0.5/W_s);         % In deg
GM          = 20*log10( 1+1/W_s );                   % In dB

% [INFO] ...
fprintf( '\t\t > PM = %2.2f deg, GM = %2.2f dB\n', PM, GM );
fprintf( '\t\t > ' );
fprintf( ACK );

%% Step 5: Define Performance Specifications

% [INFO] ...
fprintf( 'Step 5:' );
fprintf( '\tDefining performance specifications...' );

% --- Type 3
% Frequencies of interest
omega_3 = [ 0.1 0.5 1 5 10 50 ];

% Restriction
num     = [ 0.025   , 0.2   , 0.018 ];
den     = [ 0.025   , 10    , 1     ];
del_3   = tf( num, den );

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
    xlim( [-360 0] ); ylim( [-10 30] );
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
    plotbnds(bdb2);
    title('Sensitivity Reduction Bounds');
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );


%% Step 8: Intersection of QFT Bounds and Compatibility

% [INFO] ...
fprintf( 'Step 8:' );
fprintf( '\tGrouping bounds...' );

% --- Grouping bounds
bdb = grpbnds( bdb1, bdb2 );
if( PLOT )
    plotbnds(bdb); 
    title('All Bounds');
end

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

% --- Find bound intersections
ubdb = sectbnds(bdb);
if( PLOT )
    plotbnds(ubdb);
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
% --- Pole controller, G_theta(s)
G_file  = [ src 'linearInvPend_Pole_V2.shp' ];
if( isfile(G_file) )
    G = getqft( G_file );
else
    % From PID TUNER
    PID_P   = -433.469;
    PID_I   = -976.195;
    PID_D   = - 47.263;
    PID_N   =  262.366;

    % Convert to proper form
    Kp  = PID_P;
    Ti  = Kp/PID_I;
    Td  = Kp/PID_D;
    N   = PID_N;
    syms s;
    num = Kp .* sym2poly( Ti*Td*(1+1/N)*s^2 + (Ti+Td/N)*s + 1 );    % Get coefficients
    den = Ti .* sym2poly( s*( (Td/N)*s + 1 ) );                     % ...
    clear s;
    
    % Construct controller TF
    G = tf( num, den );
end

% Define a frequency array for loop shaping
wl = logspace( log10(0.01), log10(500), 2048 );
L0 = P( 1, 1, nompt );
L0.ioDelay = 0; % no delay
lpshape( wl, ubdb, L0, G );

% [INFO] ...
fprintf( ACK );

%% Step 10: Synthesize Prefitler F(s)

% [INFO] ...
fprintf( 'Step 10:' );
fprintf( '\tSynthesize F(s)...' );

syms s;
num = 1;
den = sym2poly( s/10 + 1 );
clear s;

F = tf( num, den );

pfshape( 3, wl, del_1, P, [], G, [], F );

% [INFO] ...
fprintf( ACK );

%% Step 11-13: ANALYSIS

disp(' ')
disp('chksiso(1,wl,del_1,P,R,G); %margins spec')
chksiso( 1, wl, del_1, P, [], G, [], F );
% ylim( [0 3.5] );

disp(' ')
disp('chksiso(2,wl,del_3,P,R,G); %Sensitivity reduction spec')
ind = find(wl <= 50);
chksiso( 2, wl(ind), del_3, P, [], G, [], F );
ylim( [-90 10] );


%% Check system/controller against Nyquist stability guidelines

% --- NOTE:
%   * Adding a zero corresponds to a +ve phase gain of +45deg / decade
%   * Adding a pole corresponds to a -ve phase drop of -45deg / decade
%
%   * Adding a  differentiator  shifts initial phase by +90deg
%   * Adding an integrator      shifts initial phase by -90deg
%
%   * For complex roots, phase gain/drop is +/-90deg

% Open-loop TF
T_OL = P_0*G;
[~, phi_L0] = bode( T_OL, 1e-16 );
[~, phi_Lw] = bode( T_OL, 1e+16 );
delta       = sign( phi_L0 - phi_Lw );      % +ve if Lw goes initially to the left

% Closed-loop TF
T_CL = T_OL/(1+T_OL);

% Check if Nyquist stability criterions are met
nyquistStability( tf(T_OL), false )
zpk( T_OL )

% Plot
if( PLOT )
    % Draw bode plot for further analysis
    figure();    bode( T_OL ); grid on;
    figure(); impulse( T_CL ); grid on;
end

%% Check plant against Nyquist stability guidelines

output = nyquistStability( P_0 );

if( PLOT )
    figure();  rlocus( P_0 ); grid on;
    figure(); nichols( P_0 ); grid on;
    figure(); nyquist( P_0 );
end

%% MISC. TEMPORARY OPERATIONS

clc;
% Open-loop TF
T_OL = P_0*G;
% Closed-loop TF
T_CL = T_OL/(1+T_OL);
fprintf( "\n-> G(s)\n" ); nyquistStability( tf(G), false )
fprintf( "\n-> P(s)\n" ); nyquistStability( P_0, false )
fprintf( "\n-> L(s)\n" ); nyquistStability( T_OL, false )

% Check
chksiso( 1, wl, del_1, P, [], G );

% Check impulse response
figure(); impulse( T_CL ); grid on;
