% Linearized turbine QFT control - SISO Regime 3
%   Design controller for regime 2 and regime 3 separately then combine
%   using MIMO methodology
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Aug. 20th, 2023
%
% CHANGELOG :
%   Aug. 20th, 2023
%       - Initial script
%
%   Sep. 19th, 2023
%       - Better handling of path generation for controllers
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
% PLOT = false;                               % COMMENT OUT TO PRINT FIGURES

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
src = fullfile( pathParts{1:end-2} );

% If on a UNIX machine (i.e. macOS, Ubuntu, etc...), fix path since
% strsplit() removes the leading '/' from the path.
if( isunix )
    src = [ filesep src ];
end

% Add QFT2 to path
addpath( genpath(src) );

% --- Controllers and pre-filter directories
scriptDir       = mfilename( 'fullpath' );
scriptPathParts = strsplit( scriptDir, filesep );
ctrlSrc         = fullfile( scriptPathParts{1:end-1}, 'controllerDesigns' );

% --- Directory where QFT generated controllers are stored
dirG = fullfile( ctrlSrc, 'R3', 'G' );
dirF = fullfile( ctrlSrc, 'R3', 'F' );

%% Read A, B, C, D matrices from linearized model
data_dir    = './data/';
% name_mdl    = 'SS_linearizedTurbine_MIMO_R3.mat';
name_mdl    = 'SS_linearizedTurbine_SISO_R3_inputted_GenTrq.mat';
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
                "blade120_phi"  , "blade120_omega"  , "blade120_Mact" , ...
                "blade0_phi"    , "blade0_omega"    , "blade0_Mact"   , ...
                "blade240_phi"  , "blade240_omega"  , "blade240_Mact" , ...
                "V_{wind_x}"    , "V_{wind_y}"      , "V_{wind_z}"   ];
inputNames  = [ "u_{CPC}", "u_{GenTrq}" ];
outputNames = [ "\omega_{rot}" ];
% State-space model
sys         = ss( A, B, C, D                , ...
                  'StateName' , stateNames  , ...
                  'InputName' , inputNames  , ...
                  'OutputName', outputNames );
% --- Generate TF from SS model
TF = tf( sys );

% --- Generate TF using theory
syms s;
I = eye( size(A) );
P_manual = C*(s*I - A)^-1*B + D;

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

% Variables we want to vary (Add variations)
min_A2_1  = -5.79745e-06;   max_A2_1  = 4.04086e-05 ;   grid_A2_1  = 2;
min_A2_2  = -0.0611078  ;   max_A2_2  = 0.0501627   ;   grid_A2_2  = 2;
min_A2_3  = 0.0695932   ;   max_A2_3  = 0.1002140   ;   grid_A2_3  = 2;
min_A2_5  = -0.1209990  ;   max_A2_5  = -0.0197135  ;   grid_A2_5  = 2;
min_A2_6  = 0.0716087   ;   max_A2_6  = 0.1031170   ;   grid_A2_6  = 2;
min_A2_8  = -0.1264400  ;   max_A2_8  = -0.0253774  ;   grid_A2_8  = 2;
min_A2_9  = 0.0735982   ;   max_A2_9  = 0.1059810   ;   grid_A2_9  = 2;
min_A2_11 = -0.1326090  ;   max_A2_11 = -0.03165510 ;   grid_A2_11 = 2;
min_A2_12 = 0.000950601 ;   max_A2_12 = 0.00663399  ;   grid_A2_12 = 2;


% --- Gridding
%   ***NOTE: Can grid using logspace() or linspace()
%   _g  : Gridded variable
%
A2_1_g  = linspace( (min_A2_1)    ,   (max_A2_1)  ,   grid_A2_1 );
A2_2_g  = linspace( (min_A2_2)    ,   (max_A2_2)  ,   grid_A2_2 );
A2_3_g  = linspace( (min_A2_3)    ,   (max_A2_3)  ,   grid_A2_3 );
A2_5_g  = linspace( (min_A2_5)    ,   (max_A2_5)  ,   grid_A2_5 );
A2_6_g  = linspace( (min_A2_6)    ,   (max_A2_6)  ,   grid_A2_6 );
A2_8_g  = linspace( (min_A2_8)    ,   (max_A2_8)  ,   grid_A2_8 );
A2_9_g  = linspace( (min_A2_9)    ,   (max_A2_9)  ,   grid_A2_9 );
A2_11_g = linspace( (min_A2_11)   ,   (max_A2_11) ,   grid_A2_11);
A2_12_g = linspace( (min_A2_12)   ,   (max_A2_12) ,   grid_A2_12);

% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 300 ) == SISO with 300 TFs
%
n_Plants = grid_A2_1*grid_A2_2*grid_A2_3*grid_A2_5*grid_A2_6*...
           grid_A2_8*grid_A2_9*grid_A2_11*grid_A2_12;   % Number of plants
P11 = tf( zeros(1,1,n_Plants) );                        % Pre-allocate memory
P12 = tf( zeros(1,1,n_Plants) );                        % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

NDX = 1;                                            % Plant counter
for var1 = 1:grid_A2_1                               % Loop over w
    A2_1 = A2_1_g( var1 );                            % ....
    
    for var2 = 1:grid_A2_2                           % Loop over w
        A2_2 = A2_2_g( var2 );                        % ....
        
        for var3 = 1:grid_A2_3                       % Loop over w
            A2_3 = A2_3_g( var3 );                    % ....
            
            for var4 = 1:grid_A2_5                   % Loop over w
                A2_5 = A2_5_g( var4 );                % ....
                
                for var5 = 1:grid_A2_6               % Loop over w
                    A2_6 = A2_6_g( var5 );            % ....
                    
                    for var6 = 1:grid_A2_8           % Loop over w
                        A2_8 = A2_8_g( var6 );        % ....

                        for var7 = 1:grid_A2_9       % Loop over w
                            A2_9 = A2_9_g( var7 );    % ....

                            for var8 = 1:grid_A2_11  % Loop over w
                                A2_11 = A2_11_g( var8 );        % ....

                                for var9 = 1:grid_A2_12  % Loop over w
                                    A2_12 = A2_12_g( var9 );        % ....

                                    % --- Here we create the plant TF
                                    A_g = A;    B_g = B;
                                    C_g = C;    D_g = D;
                                
                                    % Add uncertainty
                                    A_g(2, 1) = A2_1;
                                    A_g(2, 2) = A2_2;
                                    A_g(2, 3) = A2_3;
                                    A_g(2, 5) = A2_5;
                                    A_g(2, 6) = A2_6;
                                    A_g(2, 8) = A2_8;
                                    A_g(2, 9) = A2_9;
                                    A_g(2,11) = A2_11;
                                    A_g(2,12) = A2_12;
                                
                                    % --- Generate grided TF from grided SS model
                                    sys_g = ss( A_g, B_g, C_g, D_g );
                                    TF_g = tf( sys_g );
                                    P11(:, :, NDX) = TF_g(1);       % Plant TF 1,1
                                    NDX = NDX + 1;                  % Incerement counter
                                end
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
P0_11(1, 1, 1) = TF(1);                     % Nominal Transfer Function

% --- Append to the end of the gridded plants
P11( 1, 1, end+1 ) = P0_11;

% --- Cleanup plants transfer function by removing values below 1e-16
for ii = 1:length( P11 )
    [n, d] = tfdata( minreal(P11( 1, 1, ii ), 0.01) );
    n = cellfun(@(x) {x.*(abs(x) > 1e-16)}, n);
    d = cellfun(@(x) {x.*(abs(x) > 1e-16)}, d);
    P11( 1, 1, ii ) = tf(n, d);
end


% --- Define nominal plant case
nompt = length( P11 );
P0_11 = P11( 1, 1, nompt );

% --- Pick one (for now)
P0 = P0_11;
P = P11;
% P0 = P0_12;
% P = P12;

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = logspace( log10(1e-2), log10(1e1), 1024 );
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
% w = [ 1e-1 2.5e-1 5e-1 1e0 2.5e0 5e0 1e1 ];
% w = [ 5e-1 1e0 2.5e0 5e0 1e1 2.5e1 5e1 ];
w = [ 1e0 2.5e0 5e0 1e1 2.5e1 5e1 ];

% --- Plot QFT templates
if( PLOT )
    plottmpl( w, P, nompt );
    title( 'Plant Templates' )

    % --- Change legend position
    hLegend = findobj(gcf, 'Type', 'Legend');   % Get legend property
    set( hLegend, 'location', 'southeast' );    % Access and change location
    
    % --- Change plot limits
%     xmin = -405; xmax = -135; dx = 45;
%     xlim( [xmin xmax] );
%     xticks( xmin:dx:xmax )
    
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

% --------------------------------------------------
% ----      Type 1: Stability specification     ----
% --------------------------------------------------
% Frequencies of interest
% omega_1 = [ 1e-1 2.5e-1 5e-1 1e0 2.5e0 5e0 1e1 ];
% omega_1 = [ 5e-1 1e0 2.5e0 5e0 1e1 2.5e1 5e1 ];
omega_1 = [ 1e0 2.5e0 5e0 1e1 2.5e1 5e1 ];

% Restriction
W_s         = 1.66;
% W_s         = 1.46;
% W_s         = 1.08;
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

%%
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
% omega_3 = [ 5e-1 1e0 2.5e0 5e0 ];
omega_3 = [ 1e0 2.5e0 5e0 ];

% Restriction
a_d     = 5e0;
num     = [ 1/a_d   , 0 ];
den     = [ 1/a_d   , 1 ];
del_3   = tf( num, den );

% --- Plot bounds
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( del_3, min(omega_3):0.001:max(omega_3) );
    make_nice_plot();
end

%%
% --------------------------------------------------------------------
% ---- Type 4: Disturbance rejection at plant input specification ----
% --------------------------------------------------------------------
%

% Frequencies of interest
% omega_4 = [ 5e-1 1e0 2.5e0 5e0 1e1 2.5e1 5e1 ];
omega_4 = [ 1e0 2.5e0 5e0 1e1 2.5e1 5e1 ];

% Restriction
% del_4   = 0.5;
% a_U = 0.01; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.00;
a_U = 0.5; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.00;
num = [ conv([1/a_U 1], [0 1+eps_U]) ];
den = [ (1/wn)^2 (2*zeta/wn) 1 ];
% num     = [ 1/a_d   , 0 ];
% den     = [ 1/a_d   , 1 ];
del_4   = tf( num, den );

% --- Plot bounds
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( del_4, min(omega_4):0.001:max(omega_4) );
    make_nice_plot();
end

%%
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
% omega_6 = [ 1e-1 2.5e-1 5e-1 1e0 ];
% omega_6 = [ 5e-1 1e0 2.5e0 5e0 ];
omega_6 = [ 1e0 2.5e0 5e0 ];

% Restriction
% Upper bound
% a_U = 1; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.025;
a_U = 2; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.025;
num = [ conv([1/a_U 1], [0 1+eps_U]) ];
den = [ (1/wn)^2 (2*zeta/wn) 1 ];
del_6_hi = tf( num, den );
% Lower bound
% a_L = 2.5; eps_L = 0.025;
a_L = 4; eps_L = 0.025;
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
%     xlim( [-360 0] ); ylim( [-10 30] );
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

%%
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
if( PLOT )
    % [INFO] ...
    fprintf( 'Plotting bounds...' );
    
    plotbnds( bdb3 );
    title( 'Input Disturbance Rejection Bounds' );
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

%%
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

% --- Controller, G(s)
G_file  = fullfile( dirG, 'G_R3_inputted_GenTrq.shp' );
if( isfile(G_file) )
    G = getqft( G_file );
else
    syms s;
    num = (-83.33).*sym2poly( (s + 0.015) );    % Numerator
    den =           sym2poly( (s - 0.005) );    % Denominator
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

% --- Pre-filter file, F(s)
F_file  = fullfile( dirF, 'F_R3_inputted_GenTrq.fsh' );
if( isfile(F_file) )
    F = getqft( F_file );
else
    % --- Slow response filter defined below
    syms s;
    num = 0.0075;                           % Numerator
    den = sym2poly( (s + 0.0075) );         % Denominator
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

% disp(' ')
% disp('chksiso(7,wl,W3,P,R,G); %input disturbance rejection spec')
% ind = (min(omega_6) <= wl) & (wl <= max(omega_6));
% chksiso( 7, wl(ind), del_6, P, [], G, [], F );
% % ylim( [-0.1 1.3] );

%% Impulse simulations
% Some variables are manually generated, running this section as-is will
% result in errors being raised

figure();
impulse( feedback(P0*G, 1) ); grid on; hold on;
% impulse( feedback(P0*GG, 1) ); 
% impulse( feedback(P0*GGG, 1) );
title( "Impulse response" );
make_nice_plot();

figure();
step( feedback(P0*G, 1) ); grid on; hold on;
% step( feedback(P0*GG, 1) );
% step( feedback(P0*GGG, 1) );
title( "Step Response" );
make_nice_plot();
