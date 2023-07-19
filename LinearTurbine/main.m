% Linearized turbine QFT control
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Jul. 19th, 2023
%
% CHANGELOG :
%   Jul. 19th, 2023
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

%% Flags/Constants

% --- Figure counter
CNTR = 1;                                   % Figure handle counter

% --- Enable/disable plotting figures
PLOT = true;                                %#ok<NASGU> If true, plot figures!
PLOT = false;                               % COMMENT OUT TO PRINT FIGURES

% --- Enable/disable printing figures
PRNT = true;                                %#ok<NASGU>
% PRNT = false;                               % COMMENT OUT TO PRINT FIGURES

% --- [INFO] Strings
ACK = 'COMPLETED\n\n';

% --- Plot line color/style
c_line = [ 'r', 'g', 'b', 'c', 'm', 'k' ];
C_line = [ c_line, c_line, c_line, c_line ];
C_line = [ C_line, C_line, C_line, C_line ];

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

%% Read A, B, C, D matrices from linearized model
data_dir    = './data/';
name_mdl    = 'SS_linearizedTurbine.mat';
stateSpace  = load( [data_dir name_mdl ] );

% --- Get number of states
nStates = stateSpace.nx;
% --- Extract matrices from the one, big ABCD matrix
A_full = stateSpace.ABCD( 1:nStates      , 1:nStates     );
B_full = stateSpace.ABCD( 1:nStates      , nStates+1:end );
C_full = stateSpace.ABCD( nStates+1:end  , 1:nStates     );
D_full = stateSpace.ABCD( nStates+1:end  , nStates+1:end );

% --- In this case, we only care about the first 4 states
nStatesKeep = 2;
A = A_full( 1:nStatesKeep   , 1:nStatesKeep );
B = B_full( 1:nStatesKeep   , 1:end         );
C = C_full( 1:end           , 1:nStatesKeep );
D = D_full( 1:height(C)     , 1:end         );

% --- Generate state-space model
% States and inputs names
stateNames  = [ "phi" "omega" ];
inputNames  = [ "u_pitch" ];
outputNames = [ "omega" ];
% State-space model
sys         = ss( A, B, C, D                , ...
                  'StateName' , stateNames  , ...
                  'InputName' , inputNames  , ...
                  'OutputName', outputNames );

% --- Generate TF from SS model
TF = tf( sys );

%% DOUBLE CHECK
% Generate SS model using analytical approach

M_0 = 0.5   ;                   % Mass of cart                  [  kg  ]
m_0 = 0.2   ;                   % Mass of rod                   [  kg  ]
b_0 = 0.1   ;                   % Co-efficient of friction
I_0 = 0.006 ;                   % 2nd mass moment of inertia    [kg.m^2]
g_0 = 9.80665 ;                 % Gravitational acceleration    [m.s^-2]
l_0 = 0.6/2 ;                   % Half the rod length           [  m   ]

p = I_0*(M_0+m_0)+M_0*m_0*l_0^2;    % Denominator for the A and B matrices

A_theory = [ 0              1                       0               0   ;
             0  -(I_0+m_0*l_0^2)*b_0/p  (m_0^2*g_0*l_0^2)/p         0   ;
             0              0                       0               1   ;
             0  -(m_0*l_0*b_0)/p        m_0*g_0*l_0*(M_0+m_0)/p     0 ] ;
B_theory = [        0            ;
             (I_0+m_0*l_0^2)/p   ;
                    0            ;
                m_0*l_0/p       ];
C_theory = [ 1  0   0   0  ;
             0  0  -1   0 ];
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
loVal   = 0.95;             % min_ val is 95%  of nominal
hiVal   = 1.05;             % max_ val is 105% of nominal
min_M   = M_0*loVal;    max_M   = M_0*hiVal;    grid_M  = 5;
min_m   = m_0*loVal;    max_m   = m_0*hiVal;    grid_m  = 3;
min_b   = b_0*loVal;    max_b   = b_0*hiVal;    grid_b  = 5;
min_I   = I_0*loVal;    max_I   = I_0*hiVal;    grid_I  = 4;


% --- Gridding
%   ***NOTE: Can grid using logspace() or linspace()
%   _g  : Gridded variable
%
M_g = logspace( log10(min_M)    ,   log10(max_M)    ,   grid_M );
m_g = logspace( log10(min_m)    ,   log10(max_m)    ,   grid_m );
b_g = logspace( log10(min_b)    ,   log10(max_b)    ,   grid_b );
I_g = logspace( log10(min_I)    ,   log10(max_I)    ,   grid_I );


% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 300 ) == SISO with 300 TFs
%
n_Plants = grid_M*grid_m*grid_b*grid_I;             % Number of plants
P_y1 = tf( zeros(1,1,n_Plants) );                   % Pre-allocate memory
P_y2 = tf( zeros(1,1,n_Plants) );                   % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

l = l_0;                                            % Length  is nominal
g = g_0;                                            % Gravity is nominal
NDX = 1;                                            % Plant counter
for var1 = 1:grid_M                                 % Loop over M
    M = M_g( var1 );                                % ....

    for var2 = 1:grid_m                             % Loop over m
        m = m_g( var2 );                            % ....

        for var3 = 1:grid_b                         % Loop over b
            b = b_g( var3 );                        % ...

            for var4 = 1:grid_I                     % Loop over I
                I = I_g( var4 );                    % ...

                % --- Here we create the plant TF
                A_g = [ 0      1               0           0   ;
                        0  -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0   ;
                        0      0               0           1   ;
                        0  -(m*l*b)/p       m*g*l*(M+m)/p  0 ] ;
                B_g = [     0       ;
                        (I+m*l^2)/p ;
                            0       ;
                            m*l/p  ];
                C_g = [ 1 0 0 0  ;
                        0 0 -1 0 ];
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
                NDX = NDX + 1;                              % Incerement counter
            end
        end
    end
end

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
% --- Define nominal plant case
nompt = length( P );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = logspace( -4, 3 );
figure( CNTR ); CNTR = CNTR + 1;
bode( P_0, w ); grid on;
[p0, theta0] = bode( P_0, w );

make_nice_plot();

% --- Plot root locus
figure( CNTR ); CNTR = CNTR + 1;
rlocus( P_0 );
title('Root Locus of Plant (under Proportional Control)')

make_nice_plot();

%% Step 3: QFT Template

% [INFO] ...
fprintf( 'Step 3:' );
fprintf( '\tPlotting QFT templates...' );

% --- Working frequencies
% w = linspace( 1e1, 1e3, 10 );
w = [ 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 ];

% --- Plot QFT templates
plottmpl( w, P, nompt );

% --- Change legend position
hLegend = findobj(gcf, 'Type', 'Legend');   % Get legend property
set( hLegend, 'location', 'southeast' );    % Access and change location

% --- Change plot limits
xmin = -25; xmax = 10; dx = 5;
% xlim( [xmin xmax] );
% xticks( xmin:dx:xmax )
title( 'Plant Templates' )

% --- Beautify plot
make_nice_plot();

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
omega_1 = [ 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 ];
% Restriction
W_s         = 1.08;
del_1       = W_s;
PM          = 180 -2*(180/pi)*acos(0.5/W_s);         % In deg
GM          = 20*log10( 1+1/W_s );                   % In dB

% [INFO] ...
fprintf( '\t\t > PM = %2.2f deg, GM = %2.2f dB', PM, GM );
fprintf( '\t\t > ' );
fprintf( ACK );

%% Step 5: Define Performance Specifications

% [INFO] ...
fprintf( 'Step 5:' );
fprintf( '\tDefining performance specifications...' );

% --- Type 3
% Frequencies of interest
omega_3 = [ 1e-2 1e-1 1e0 1e1 ];

% Restriction
num     = [ 0.025   , 0.2   , 0.018 ];
den     = [ 0.025   , 10    , 1     ];
del_3   = tf( num, den );


% --- Type 6
% Frequencies of interest
omega_6 = [ 1e-2 1e-1 1e0 1e1 ];

% Restriction
% Upper bound
a_U = 0.1; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.05;
num = [ 21/2 21/20 ];
den = [ (1/wn)*(1/wn) (2*zeta/wn) 1 ];
del_6_hi = tf( num, den );
% Lower bound
a_L = 0.25; eps_L = 0.0;
num = 1-eps_L;
den = [ 16 8 1 ];
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

% [INFO] ...
fprintf( 'Plotting bounds...' );

% --- Plot bounds
plotbnds( bdb1 );
title( 'Robust Stability Bounds' );
xlim( [-360 0] ); ylim( [-10 30] );
make_nice_plot();

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

% [INFO] ...
fprintf( 'Plotting bounds...' );

% --- Plot bounds
plotbnds(bdb2);
title('Sensitivity Reduction Bounds');
make_nice_plot();

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

% [INFO] ...
fprintf( 'Plotting bounds...' );

% --- Plot bounds
plotbnds(bdb7);
title('Robust Tracking Bounds');
make_nice_plot();

% [INFO] ...
fprintf( ACK );


%% Step 8: Intersection of QFT Bounds and Compatibility

% [INFO] ...
fprintf( 'Step 8:' );
fprintf( '\tGrouping bounds...' );

% --- Grouping bounds
bdb = grpbnds( bdb1, bdb2 );
% bdb = grpbnds( bdb1, bdb2, bdb7 );
plotbnds(bdb); 
title('All Bounds');

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

% --- Find bound intersections
ubdb = sectbnds(bdb);
plotbnds(ubdb);
title('Intersection of Bounds');

% [INFO] ...
fprintf( ACK );

%% Step 9: Synthesize Feedback Controller G(s)

% [INFO] ...
fprintf( 'Step 9:' );
fprintf( '\tSynthesize G(s)...' );

% --- Design TF of G(s)
syms s;
num = sym2poly( 250*(s/1 + 1)*(s/1000 + 1)  );      % Get coefficients
den = sym2poly( s*(s/5 + 1)*(s/1000 + 1));          % ...
clear s;
% syms s;
% num = sym2poly( 1000*(s/55 + 1)*(s/1700 + 1)  );  % Get coefficients
% den = sym2poly( s*(s/25000 + 1)               );  % ...
% clear s;

% Construct controller TF
G = tf( num, den );

% Define a frequency array for loop shaping
wl = linspace( 1, 500000 );
L0 = P( 1, 1, nompt );
L0.ioDelay = 0; % no delay
lpshape( wl, ubdb, L0, G );

% [INFO] ...
fprintf( ACK );


%% Step 10: Synthesize Prefitler F(s)
%{
% [INFO] ...
fprintf( 'Step 10:' );
fprintf( '\tSynthesize F(s)...' );

syms s;
num = 1;
den = sym2poly( s/550 + 1 );
clear s;

F = tf( num, den );

pfshape( 7, wl, del_6, P, [], G, [], F );

% [INFO] ...
fprintf( ACK );

%% Step 11-13: ANALYSIS

disp(' ')
disp('chksiso(1,wl,del_1,P,R,G); %margins spec')
chksiso( 1, wl, del_1, P, [], G );
ylim( [0 3.5] );

disp(' ')
disp('chksiso(2,wl,del_3,P,R,G); %Sensitivity reduction spec')
ind = find(wl <= 1000);
chksiso( 2, wl(ind), del_3, P, [], G );
ylim( [-90 10] );

disp(' ')
disp('chksiso(7,wl,W3,P,R,G); %input disturbance rejection spec')
drawnow
% chksiso(7,wl(ind),W7,P,[],G);
chksiso( 7, wl, del_6, P, [], G, [], F );
% ylim( [-0.1 1.3] );
%}