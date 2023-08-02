% Linearized turbine QFT control
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Jul. 19th, 2023
%
% CHANGELOG :
%   Jul. 19th, 2023
%       - Initial script
%
%   Aug.  2nd, 2023
%       - Added generator torque for regime 3 operations
%       - Added blade pitching actuator dynamics
%       - Added wind disturbance
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
% name_mdl    = 'SS_linearizedTurbine_16p46rpm.mat';
name_mdl    = 'SS_linearizedTurbine_Rated_w_GenTrq_w_ActuatorDynamics.mat';
stateSpace  = load( [data_dir name_mdl ] );

% --- Get number of states
nStates = stateSpace.nx;
% --- Extract matrices from the one, big ABCD matrix
A_full = stateSpace.ABCD( 1:nStates      , 1:nStates     );
B_full = stateSpace.ABCD( 1:nStates      , nStates+1:end );
C_full = stateSpace.ABCD( nStates+1:end  , 1:nStates     );
D_full = stateSpace.ABCD( nStates+1:end  , nStates+1:end );

% --- In this case, we only care about the first "nStatesKeep" states
nStatesKeep = nStates;
A = A_full( 1:nStatesKeep   , 1:nStatesKeep );
B = B_full( 1:nStatesKeep   , 1:end         );
C = C_full( 1:end           , 1:nStatesKeep );
D = D_full( 1:height(C)     , 1:end         );

% --- Generate state-space model
% States and inputs names
stateNames  = [ "phi"           , "omega"           , ...
                "blade120_phi"  , "blade120_omega"  , "blade120_Mact", ...
                "blade0_phi"    , "blade0_omega"    , "blade0_Mact"  , ...
                "blade240_phi"  , "blade240_omega"  , "blade2400_Mact"];
inputNames  = [ "u_{pitch}" ];
outputNames = [ "omega" ];
% State-space model
sys         = ss( A, B, C, D                , ...
                  'StateName' , stateNames  , ...
                  'InputName' , inputNames  , ...
                  'OutputName', outputNames );
% --- Generate TF from SS model
TF = tf( sys );

%% __DEPRECATED__Actuator Dynamics, M_act(s)

% --- The blade pitching actuator has a rate of 10 deg/s <==> 0.1745 rad/s
%
pitch_rate_deg = 10;                        % Pitching rate in [deg/s]
pitch_rate_rad = pitch_rate_deg*pi/180;     % Pitching rate in [rad/s]
pitch_rate_hz  = pitch_rate_rad/(2*pi);     % Pitching rate in [  Hz ]
M_act = tf( 1, [pitch_rate_rad 1] );

%% Step 1: Plant Modeling & Uncertainty

%   For this specific case, disturbances in wind speed and rotor angular
% velocity cause changes in the matrix elements listed below. Therefore, we
% add the uncertainties into those elements
%

% --- Plant parameters
%   min_    : Minimum value
%   max_    : Maximum value
%   grid_   : Gridding
%
loVal   = 0.90;             % min_ val is 90%  of nominal
hiVal   = 1.10;             % max_ val is 110% of nominal
% Variables we want to vary
A21     = A(2, 1);  A22     = A(2, 2);
A23     = A(2, 3);
A25     = A(2, 5);  A26     = A(2, 6);
A28     = A(2, 8);  A29     = A(2, 9);
A211    = A(2,11);
% Add variations
min_A21 = A21*loVal;    max_A21 = A21*hiVal;    grid_A21 = 2;
min_A22 = A22*loVal;    max_A22 = A22*hiVal;    grid_A22 = 2;
min_A23 = A23*loVal;    max_A23 = A23*hiVal;    grid_A23 = 2;
min_A25 = A25*loVal;    max_A25 = A25*hiVal;    grid_A25 = 2;
min_A26 = A26*loVal;    max_A26 = A26*hiVal;    grid_A26 = 2;
min_A28 = A28*loVal;    max_A28 = A28*hiVal;    grid_A28 = 2;
min_A29 = A29*loVal;    max_A29 = A29*hiVal;    grid_A29 = 2;
min_A211= A211*loVal;   max_A211= A211*hiVal;   grid_A211= 2;


% --- Gridding
%   ***NOTE: Can grid using logspace() or linspace()
%   _g  : Gridded variable
%
A21_g = linspace( (min_A21)    ,   (max_A21)  ,   grid_A21 );
A22_g = linspace( (min_A22)    ,   (max_A22)  ,   grid_A22 );
A23_g = linspace( (min_A23)    ,   (max_A23)  ,   grid_A23 );
A25_g = linspace( (min_A25)    ,   (max_A25)  ,   grid_A25 );
A26_g = linspace( (min_A26)    ,   (max_A26)  ,   grid_A26 );
A28_g = linspace( (min_A28)    ,   (max_A28)  ,   grid_A28 );
A29_g = linspace( (min_A29)    ,   (max_A29)  ,   grid_A29 );
A211_g= linspace( (min_A211)   ,   (max_A211) ,   grid_A211);
% A21_g = logspace( log10(min_A21)    ,   log10(max_A21)  ,   grid_A21 );
% A22_g = logspace( log10(min_A22)    ,   log10(max_A22)  ,   grid_A22 );
% A23_g = logspace( log10(min_A23)    ,   log10(max_A23)  ,   grid_A23 );

% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 300 ) == SISO with 300 TFs
%
n_Plants = grid_A21*grid_A22*grid_A23*grid_A25*...
           grid_A26*grid_A28*grid_A29*grid_A211;    % Number of plants
P = tf( zeros(1,1,n_Plants) );                      % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

NDX = 1;                                            % Plant counter
for var1 = 1:grid_A21                               % Loop over w
    A21 = A21_g( var1 );                            % ....
    
    for var2 = 1:grid_A22                           % Loop over w
        A22 = A22_g( var2 );                        % ....
        
        for var3 = 1:grid_A23                       % Loop over w
            A23 = A23_g( var3 );                    % ....
            
            for var4 = 1:grid_A25                   % Loop over w
                A25 = A25_g( var4 );                % ....
                
                for var5 = 1:grid_A26               % Loop over w
                    A26 = A26_g( var5 );            % ....
                    
                    for var6 = 1:grid_A28           % Loop over w
                        A28 = A28_g( var6 );        % ....

                        for var7 = 1:grid_A29       % Loop over w
                            A29 = A29_g( var7 );    % ....

                            for var8 = 1:grid_A211  % Loop over w
                                A211 = A211_g( var8 );        % ....

                                % --- Here we create the plant TF
                                A_g = A;    B_g = B;
                                C_g = C;    D_g = D;
                            
                                % Add uncertainty
                                A_g(2, 1) = A21;
                                A_g(2, 2) = A22;
                                A_g(2, 3) = A23;
                                A_g(2, 5) = A25;
                                A_g(2, 6) = A26;
                                A_g(2, 8) = A28;
                                A_g(2, 9) = A29;
                                A_g(2,11) = A211;
                            
                                % --- Generate grided TF from grided SS model
                                sys_g = ss( A_g, B_g, C_g, D_g );
                                TF_g = tf( sys_g );
                                P(:, :, NDX) = TF_g(1);         % Transfer Function
                                NDX = NDX + 1;                  % Incerement counter
                            end
                        end
                    end
                end
            end
        end
    end
end 

% [INFO] ...
fprintf( ACK );

%% Step 2: The Nominal Plant

% [INFO] ...
fprintf( 'Step 2:' );
fprintf( '\tComputing nominal plant...' );

% --- Generate nominal plant TF
%   Any one of the models above can be used as the nominal plant.
%   We just happened to chose this one.
%
P0(1, 1, 1) = TF;                       % Nominal Transfer Function

% --- Append to the end of the gridded plants
P( 1, 1, end+1 ) = P0;

% --- Cleanup plants transfer function by removing values below 1e-16
for ii = 1:length( P )
    [n, d] = tfdata( minreal(P( 1, 1, ii ), 0.01) );
    n = cellfun(@(x) {x.*(abs(x) > 1e-16)}, n);
    d = cellfun(@(x) {x.*(abs(x) > 1e-16)}, d);
    P( 1, 1, ii ) = tf(n, d);
end

% --- Define nominal plant case
nompt = length( P );
P0 = P( 1, 1, nompt );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = logspace( log10(1e-4), log10(1e2), 1024 );
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( P0, w ); grid on;
    make_nice_plot();
end
[p0, theta0] = bode( P0, w );

% --- Plot root locus
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    rlocus( P0 );
    title('Root Locus of Plant')
    make_nice_plot();
end

%% Step 3: QFT Template

% [INFO] ...
fprintf( 'Step 3:' );
fprintf( '\tPlotting QFT templates...' );

% --- Working frequencies
% w = linspace( 5e-3, 1e1, 8 );
w = [ 5e-3 1e-2 5e-2 1e-1 5e-1 1e0 5e0 1e1 ];

% --- Plot QFT templates
plottmpl( w, P, nompt );

% --- Change legend position
hLegend = findobj(gcf, 'Type', 'Legend');   % Get legend property
set( hLegend, 'location', 'southeast' );    % Access and change location

% --- Change plot limits
% xmin = -315; xmax = -135; dx = 45;
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

% --------------------------------------------------
% ----      Type 1: Stability specification     ----
% --------------------------------------------------
% Frequencies of interest
% omega_1 = [ 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 ];
omega_1 = [ 1e-2 5e-2 1e-1 5e-1 1e0 5e0 1e1 ];
% Restriction
% W_s         = 1.46;
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

% -----------------------------------------------------------------------
% -- Type 3: Sensitivity or output disturbance rejection specification --
% -----------------------------------------------------------------------
%
% Typically use the following form:
%   
%   --> del_3(s) = s/(s + a_d)
%
%   By selecting just one parameter, the pole a_d, we can achieve different
% levels of disturbance rejection. The higher the parameter a_d, the more
% significant the attenuation of the effect of the disturbance.
%
%   Typical choice for a_d is s.t. ( a_d >= max(omega_3) )
%

% Frequencies of interest
omega_3 = [ 1e-2 5e-2 1e-1 ];

% Restriction
a_d     = 1e-1;
num     = [ 1/a_d   , 0 ];
den     = [ 1/a_d   , 1 ];
% num     = [ 0.025   , 0.2   , 0.018 ];
% den     = [ 0.025   , 10    , 1     ];
del_3   = tf( num, den );

% --- Plot bounds
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( del_3, min(omega_3):0.01:max(omega_3) );
    make_nice_plot();
end


% --------------------------------------------------------------------
% ---- Type 4: Disturbance rejection at plant input specification ----
% --------------------------------------------------------------------
%

% Frequencies of interest
omega_4 = [ 1e-2 5e-2 1e-1 5e-1 1e0 5e0 1e1 ];

% Restriction
del_4   = 2;
% a_U = 0.10; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.01;
% num = [ conv([1/a_U 1], [0 1+eps_U]) ];
% den = [ (1/wn)^2 (2*zeta/wn) 1 ];
% num     = [ 1/a_d   , 0 ];
% den     = [ 1/a_d   , 1 ];
% del_4   = tf( num, den );

% --- Plot bounds
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( del_4, min(omega_4):0.01:max(omega_4) );
    make_nice_plot();
end

% --------------------------------------------------
% ---- Type 6: Reference tracking specification ----
% --------------------------------------------------
%
% A practical selection is:
%
%                         (1-eps_L)
%   --> del_6_lo(s) = -----------------
%                       (s/a_L + 1)^2
%       With
%               0 <= eps_L
%
%
%                           (s/a_U + 1)*(1+eps_U)
%   --> del_6_hi(s) = ----------------------------------
%                       ((s/wn)^2 + (2*zeta*s/wn) + 1)
%       With
%               0 <= eps_U; zeta = 0.8; wn = 1.25*a_U/zeta
%
%   Normally, we do not ask the system to follow a high-frequency
% reference. In this way, we reduce high-frequency activity of the
% actuators and then avoid potential mechanical fatigue problems.
%

% Frequencies of interest
omega_6 = [ 1e-2 5e-2 1e-1 1e0 ];

% Restriction
% Upper bound
a_U = 0.25; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.01;
num = [ conv([1/a_U 1], [0 1+eps_U]) ];
den = [ (1/wn)^2 (2*zeta/wn) 1 ];
del_6_hi = tf( num, den );
% Lower bound
a_L = 0.50; eps_L = 0.01;
num = 1-eps_L;
den = [ conv([1/a_L 1], [1/a_L 1]) ];
del_6_lo = tf( num, den );
% Tracking weight
del_6 = [ del_6_hi  ;
          del_6_lo ];

% --- Plot bounds
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    step( del_6(1) );   hold on ;  grid on;
    step( del_6(2) );   hold off;
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

%% Step 6: Calculate Staibility QFT Bounds

% --- Type refers to Dr. Garcia's book
% --- ptype referes to the QFT toolbox
%
%   - Type 1: Stability specification
%       > Corresponds to sisobnds( ptype=1, ... )
%   - Type 3    //
%       > Corresponds to sisobnds( ptype=2, ... )
%   - Type 4    //
%       > Corresponds to sisobnds( ptype=3, ... )
%   - Type 6    //
%       > Corresponds to sisobnds( ptype=7, ... )
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

% --- Plot bounds
if( PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );

    plotbnds( bdb1 );
    title( 'Robust Stability Bounds' );
    xlim( [-360 0] ); ylim( [-10 30] );
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

%% Step 7: Calculate Performance QFT Bounds

% -----------------------------------------------------------------------
% -- Type 3: Sensitivity or output disturbance rejection specification --
% -----------------------------------------------------------------------
spec = 2;

% [INFO] ...
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb2 = sisobnds( spec, omega_3, del_3, P, [], nompt );

% --- Plot bounds
if( PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );
    
    plotbnds( bdb2 );
    title( 'Sensitivity Reduction Bounds' );
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

% --------------------------------------------------------------------
% ---- Type 4: Disturbance rejection at plant input specification ----
% --------------------------------------------------------------------
spec = 3;

% [INFO] ...
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb3 = sisobnds( spec, omega_4, del_4, P, [], nompt );

% --- Plot bounds
if( ~PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );
    
    plotbnds( bdb3 );
    title( 'Input Disturbance Rejection Bounds' );
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

% --- Plot bounds
if( PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );
    
    plotbnds( bdb7 );
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
bdb = grpbnds( bdb1, bdb2, bdb3, bdb7 );
% --- Plot bounds
if( PLOT )
    plotbnds( bdb ); 
    title( 'All Bounds' );
end

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

% --- Find bound intersections
ubdb = sectbnds(bdb);
% --- Plot bounds
if( PLOT )
    plotbnds( ubdb );
    title( 'Intersection of Bounds' );
end

% [INFO] ...
fprintf( ACK );

%% Step 9: Synthesize Feedback Controller G(s)

% [INFO] ...
fprintf( 'Step 9:' );
fprintf( '\tSynthesize G(s)...' );

% --- Directory where QFT generated controllers are stored
src = './controllerDesigns/';

% --- Controller, G(s)
% G_file  = [ src 'G.shp' ];
G_file  = [ src 'G_w_GenTrq_w_ActuatorDynamics_w_UpdatedSpecs_w_WindDisturbance.shp' ];
if( isfile(G_file) )
    G = getqft( G_file );
else
    syms s;
    num = (-7.5) .* sym2poly( (s/0.1 + 1) );    % Numerator
    den =           sym2poly( (s/8   + 1) );    % Denominator
    clear s;
    
    % Construct controller TF
    G = tf( num, den );                         % Eq.(CS3.25)
end

% Define a frequency array for loop shaping
wl = logspace( log10(min(w)*1e-1), log10(max(w)*1e1), 1024 );
L0 = P( 1, 1, nompt );
L0.ioDelay = 0; % no delay
lpshape( wl, ubdb, L0, G );

% --- Store as SS in case we want to use a SS representation in Simulink
[A_G, B_G, C_G, D_G] = tf2ss( cell2mat(tf(G).num), cell2mat(tf(G).den) );

% [INFO] ...
fprintf( ACK );


%% Step 10: Synthesize Prefitler F(s)

% [INFO] ...
fprintf( 'Step 10:' );
fprintf( '\tSynthesize F(s)...' );

% --- Directory where QFT generated controllers are stored
src = './controllerDesigns/';
% --- Pre-filter file, F(s)
F_file  = [ src 'F_w_GenTrq_w_ActuatorDynamics_w_UpdatedSpecs_w_WindDisturbance.fsh' ];
if( isfile(F_file) )
    F = getqft( F_file );
else
    syms s;
    num = 1.05;                                     % Numerator
    den = sym2poly( (s/0.225 + 1)*(s/2 + 1) );      % Denominator
    clear s;
    
    % Construct controller TF
    F = tf( num, den );
end

pfshape( 7, min(omega_6):0.001:max(omega_6), del_6, P, [], G, [], F );

% --- Store as SS in case we want to use a SS representation in Simulink
[A_F, B_F, C_F, D_F] = tf2ss( cell2mat(tf(F).num), cell2mat(tf(F).den) );

% [INFO] ...
fprintf( ACK );

%% Step 11-13: ANALYSIS

disp(' ')
disp('chksiso(1,wl,del_1,P,R,G); %margins spec')
ind = (min(omega_1) <= wl) & (wl <= max(omega_1));
chksiso( 1, wl(ind), del_1, P, [], G );
% ylim( [0 3.5] );

disp(' ')
disp('chksiso(2,wl,del_3,P,R,G); %Sensitivity reduction spec')
ind = (min(omega_3) <= wl) & (wl <= max(omega_3));
chksiso( 2, wl(ind), del_3, P, [], G );
% ylim( [-90 10] );

disp(' ')
disp('chksiso(3,wl,del_4,P,R,G); %Disturbance at input reduction spec')
ind = (min(omega_4) <= wl) & (wl <= max(omega_4));
chksiso( 3, wl(ind), del_4, P, [], G );
% ylim( [-90 10] );

disp(' ')
disp('chksiso(7,wl,W3,P,R,G); %input disturbance rejection spec')
ind = (min(omega_6) <= wl) & (wl <= max(omega_6));
chksiso( 7, wl(ind), del_6, P, [], G, [], F );
% ylim( [-0.1 1.3] );
