% Waste Water Treatment Plant (WWTP) QFT control design
%   Based on Dr. Mario Garcia-Sanz's book
%       "Robust Control eEngineering: Practical QFT Solutions"
%       Case Study 3 - Pg. 343
%
%   2x2 MIMO system
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Jul.  6th, 2023
%
% CHANGELOG :
%   Jul.  6th, 2023
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
% PLOT = false;                               % COMMENT OUT TO PRINT FIGURES

% --- Enable/disable printing figures
PRNT = true;                                %#ok<NASGU>
% PRNT = false;                               % COMMENT OUT TO PRINT FIGURES

% --- [INFO] Strings
ACK = 'COMPLETED\n\n';

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

%% Step 1: Plant Modeling & Uncertainty

% --- Plant parameters
%   min_    : Minimum value
%   max_    : Maximum value
%   grid_   : Gridding
%
min_k11 = -0.045 ;  max_k11 = -0.035 ;  grid_k11  = 3;
min_a11 = 6.25e-5;  max_a11 = 8.25e-5;  grid_a11  = 3;
min_k22 = -2.2e-5;  max_k22 = -1.8e-5;  grid_k22  = 3;
min_a22 = 1.57e-4;  max_a22 = 1.77e-4;  grid_a22  = 3;

% --- Gridding
%   ***NOTE: Can grid using logspace() or linspace()
%   _g  : Gridded variable
%
k11_g = linspace( min_k11, max_k11, grid_k11 );
a11_g = linspace( min_a11, max_a11, grid_a11 );
k22_g = linspace( min_k22, max_k22, grid_k22 );
a22_g = linspace( min_a22, max_a22, grid_a22 );
% k11_g = logspace( log10(min_k11), log10(max_k11), grid_k11 );
% a11_g = logspace( log10(min_a11), log10(max_a11), grid_a11 );
% k22_g = logspace( log10(min_k22), log10(max_k22), grid_k22 );
% a22_g = logspace( log10(min_a22), log10(max_a22), grid_a22 );

% --- Constant parameters
k12     = -6.239e-6;
z12_1   =  7.534e-4;
z12_2   = -3.170e-5;
a12     =  8.040e-5;
wn12    =  4.580e-4;
zeta12  =  0.849300;
k21     =  0.046400;
a21     =  1.008e-4;

% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 300 ) == SISO with 300 TFs
%
n_Plants = grid_k11*grid_a11*grid_k22*grid_a22;     % Number of plants
p11 = tf( zeros(1,1,n_Plants) );                    % Pre-allocate memory
p12 = tf( zeros(1,1,n_Plants) );                    % Pre-allocate memory
p21 = tf( zeros(1,1,n_Plants) );                    % Pre-allocate memory
p22 = tf( zeros(1,1,n_Plants) );                    % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

NDX = 1;                                            % Plant counter
for var1 = 1:grid_k11                               % Loop over k11
    k11 = k11_g( var1 );                            % ....
    
    for var2 = 1:grid_a11                           % Loop over a11
        a11 = a11_g( var2 );                        % ....

        for var3 = 1:grid_k22                       % Loop over k22
            k22 = k22_g( var3 );                    % ....

            for var4 = 1:grid_a22                   % Loop over a22
                a22 = a22_g( var4 );                % ....
                
                % --- Here we create the plant TF
                p11(:, :, NDX) = tf( k11        , ...
                                    [1/a11, 1] );
                p12(:, :, NDX) = tf( k12*conv([1/z12_1, 1], [1/z12_2, 1]), ...
                                     conv([1/a12  , 1], [1/wn12^2, (2*zeta12/wn12)^2, 1]) );
                p21(:, :, NDX) = tf( k21        , ...
                                    [1/a21, 1] );
                p22(:, :, NDX) = tf( k22        , ...
                                    [1/a22, 1] );
                % % --- Place them all in one big matrix
                % Pw(:, :, NDX) = [ p11, p12 ;
                %                   p21, p22 ];
                NDX = NDX + 1;                      % Increment counter
            end
        end
    end
end

% --- Place them all in one big matrix
% ***NOTE:
%   Pw( 1, 1,  1, : ) ==> p11 / 1st  variation (i.e. p11(:,:,1))
%   Pw( 2, 1, 65, : ) ==> p21 / 65th variation (i.e. p21(:,:,65))
%
Pw = [ p11, p12 ;
       p21, p22 ];


% --- EXTRA STEP: modify Pw(s) as per the problem requirement
% Add low-ass filter
f_LP      = tf( 1, [1/(1.1e-4)^2, 2/(1.1e-4), 1 ] );
% Generate modified plants
P = Pw.*f_LP;

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
k11_0 = mean( [min_k11, max_k11] );
a11_0 = mean( [min_a11, max_a11] );
k22_0 = mean( [min_k22, max_k22] );
a22_0 = mean( [min_a22, max_a22] );

p11_0 = tf( k11_0, [1/a11_0, 1] );
p12_0 = tf( k12*conv([1/z12_1, 1], [1/z12_2, 1]), ...
            conv([1/a12  , 1], [1/wn12^2, (2*zeta12/wn12)^2, 1]) );
p21_0 = tf( k21, [1/a21, 1] );
p22_0 = tf( k22_0, [1/a22_0, 1] );

% Nominal plant TF
Pw_0 = [ p11_0, p12_0;
         p21_0, p22_0];
% Modified nominal plant
P_0 = Pw_0.*f_LP;

% --- Append to the end of the gridded plants
P( :, :, :, end+1 ) = P_0;

% --- Define nominal plant case
nompt = length( P );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = logspace( -7, -2, 1024 );
[p0, theta0] = bode( P_0, w );

if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( P_0, '-', P_0, '.r', w(1:32:end) ); grid on;
    make_nice_plot();
end

% --- Plot root locus
% if( PLOT )
%     figure( CNTR ); CNTR = CNTR + 1;
%     rlocus( P_0 );
%     title('Root Locus of Plant (under Proportional Control)')
%     make_nice_plot();
% end

%% Step 3: QFT Template

% [INFO] ...
fprintf( 'Step 3:' );
fprintf( '\tPlotting QFT templates...' );

% --- Working frequencies
w =  [ 1e-7 5e-7, 1e-6 5e-6, 1e-5 2e-5 5e-5, 1e-4 2e-4 5e-4, 1e-3 2e-4 5e-3 ];

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
G_file  = [ src 'linearInvPend_Pole_Simplified_V2.shp' ];
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
% 
% % [INFO] ...
% fprintf( 'Step 10:' );
% fprintf( '\tSynthesize F(s)...' );
% 
% syms s;
% num = 1;
% den = sym2poly( s/10 + 1 );
% clear s;
% 
% F = tf( num, den );
% 
% pfshape( 1, wl, del_1, P, [], G, [], F );
% 
% % [INFO] ...
% fprintf( ACK );

%% Step 11-13: ANALYSIS

disp(' ')
disp('chksiso(1,wl,del_1,P,R,G); %margins spec')
chksiso( 1, wl, del_1, P, [], G );
% ylim( [0 3.5] );

disp(' ')
disp('chksiso(2,wl,del_3,P,R,G); %Sensitivity reduction spec')
ind = find(wl <= 50);
chksiso( 2, wl(ind), del_3, P, [], G );
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

% Check SISO for sensitivity reduction
% ind = find(wl <= 50);
% chksiso( 2, wl(ind) , del_3, P, [], G );
chksiso( 1, wl      , del_1, P, [], G );
ylim( [-90 10] );

% Check impulse response
figure(); impulse( T_CL ); grid on;