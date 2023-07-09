% Waste Water Treatment Plant (WWTP) QFT control design
%   Based on Dr. Mario Garcia-Sanz's book
%       "Robust Control eEngineering: Practical QFT Solutions"
%       Case Study 3 - Pg. 343
%
%   2x2 MIMO system
%       Independent SISO Controller Design
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
PLOT = false;                               % COMMENT OUT TO PRINT FIGURES

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
P0 = Pw_0.*f_LP;

% --- Append to the end of the gridded plants
P( :, :, end+1, : ) = P0;

% --- Define nominal plant case
nompt = length( P );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = logspace( -7, -2, 1024 );
[p0, theta0] = bode( P0, w );

if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( P0, '-', P0, '.r', w(1:32:end) ); grid on;
    make_nice_plot();
end

%% Step 3: QFT Template

% [INFO] ...
fprintf( 'Step 3:' );
fprintf( '\tPlotting QFT templates...' );

% --- Working frequencies
w =  [ 1e-7 5e-7, 1e-6 5e-6, 1e-5 2e-5 5e-5, 1e-4 2e-4 5e-4, 1e-3 2e-4 5e-3 ];

if( PLOT )
    % --- Plot QFT templates
    for ROW = 1:width(P)
        for COL = 1:width(P)
            plottmpl( w, P(ROW, COL, :, :), nompt );
    
            % --- Change legend position
            hLegend = findobj( gcf, 'Type', 'Legend' ); % Get legend property
            set( hLegend, 'location', 'southeast' );    % Access and change location
            
            % --- Change plot limits
            if( ROW == 2 && COL == 1)
                xmin = -270; xmax = 45; dx = 45;
                xlim( [xmin xmax] );
                xticks( xmin:dx:xmax )
            end

            txt = ['Plant Templates for p' num2str(ROW) num2str(COL) '(s)' ];
            title( txt )
            
            % --- Beautify plot
            make_nice_plot();
        end
    end
end

% [INFO] ...
fprintf( ACK );

%% Step 3.5: RGA of the nominal plant, P0(s=0) and P0(s=inf)

% --- Relative Gain Array (RGA) analysis
%

% Recall, the RGA matrix, Lambda, is defined as
%
%   Λ_0   = P( s=0 ) .* (P( s=0 )^-1)^T
%               AND
%   Λ_inf = P(s=inf) .* (P(s=inf)^-1)^T
%

% --- RGA for s=jw=0
P0_0 = freqresp( P0, 0 );
Lambda_0 = P0_0 .* inv(P0_0).';

% --- RGA for s=jw=inf
P0_inf = freqresp( P0, 1e16 );
Lambda_inf = P0_inf .* inv(P0_inf).';

% --- Determine pairing
% Recall, the column element closest to 1 corresponds to the row pairing
%
%   Example:
%                                        _   u1      u2      u3   _
%                                       |  0.3180  0.0195  0.6630  | y1
%       Λ_0 = P(s=0) .* (P(s=0)^-1)^T = |  0.6820  0.0091  0.3090  | y2
%                                       |_    0    0.9710  0.0287 _| y3
%
%   According to RGA matrix, pairing is:
%       ( u1, y2 ) --- ( u2, y3 ) --- ( u3, y1 )
%

fprintf( "Control - Output pairing:\n" );   % [INFO] ...
VAL = 1;                                    % Value we want to be close to
for COL = 1:width(Lambda_0)
    Lambda_COL = Lambda_0(:, COL);          % Extract column
    % Get the index of the element closest to VAL (=1)
    [minValue, NDX] = min( abs(Lambda_COL-VAL) );
    closestValue = Lambda_COL( NDX );

    fprintf( "\t> ( u%i, y%i )\n", COL, NDX );
end

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
omega_1 = [ 1e-7 5e-7, 1e-6  5e-6, 1e-5 2e-5 5e-5, ...
            1e-4 2e-4  5e-4, 1e-3  2e-4 5e-3 ];
% Restriction (for p_ii, i=1,2)
W_s         = 1.66;
del_1       = W_s;
PM          = 180 - 2*(180/pi)*acos(0.5/W_s);       % In deg
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
%   Sensitivity or disturbances at plant output specification
%

% Frequencies of interest
omega_3 = [ 1e-7 5e-7, 1e-6  5e-6, 1e-5 2e-5 ];

% Restriction
a_d     = 2e-5;
num     = [ 1/a_d , 0 ];
den     = [ 1/a_d , 1 ];
del_3   = tf( num, den );

if( PLOT )
    % --- PLOT bode response of del_3(s)
    figure( CNTR ); CNTR = CNTR + 1;

    w_del_3 = logspace( log10(w(1)), log10(w(end)));
    [mag, ~] = bode( del_3, w_del_3 );
    mag_dB = db( squeeze(mag) );

    semilogx( w_del_3, mag_dB ); grid on;
end

% --- Type 6
%   Reference tracking specification
%

% Frequencies of interest
omega_6 = [ 1e-7 5e-7, 1e-6  5e-6, 1e-5 2e-5 5e-5, 1e-4 ];

% Restriction
% Upper bound
a_U = 2e-5; zeta = 0.8; wn = 1.25*a_U/zeta;
num = [ 1/a_U , 1 ];
den = [ 1/wn^2, 2*zeta/wn, 1 ];
del_6_U = tf( num, den );
% Lower bound
syms s;
a_L = 2e-5;
num = 1;
den = sym2poly( (s/a_L + 1)^2 );
del_6_L = tf( num, den );
clear s;
% Tracking weight
del_6 = [ del_6_U  ;
          del_6_L ];

if( PLOT )
    % --- PLOT step response of del_6_U(s) and del_6_L(s)
    figure( CNTR ); CNTR = CNTR + 1;
    stepplot( del_6_U, del_6_L ); grid on;
end

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
for i=1:width(P)
    p_ii = P( i, i, :, : );
    bdb1(:, :, i) = sisobnds( spec, omega_1, del_1, p_ii, [], nompt );
    % R = 0; bdb1 = sisobnds( spec, omega_1, del_1, P, R, nompt );
    
    if( PLOT )
        % [INFO] ...
        fprintf( 'Plotting bounds...' );
        
        % --- Plot bounds
        plotbnds( bdb1(:, :, i) );
        txt = ['Robust Stability Bounds for p' num2str(i) num2str(i) '(s)' ];
        title( txt );
        % xlim( [-360 0] ); ylim( [-10 30] );
        make_nice_plot();
    end
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
for i=1:width(P)
    p_ii = P( i, i, :, : );
    bdb2(:, :, i) = sisobnds( spec, omega_3, del_3, p_ii, [], nompt );
    
    if( PLOT )
        % [INFO] ...
        fprintf( 'Plotting bounds...' );
        
        % --- Plot bounds
        plotbnds( bdb2(:, :, i) );
        txt = ['Sensitivity Reduction Bounds for p' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot();
    end
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
for i=1:width(P)
    p_ii = P( i, i, :, : );
    bdb7(:, :, i) = sisobnds( spec, omega_6, del_6, p_ii, [], nompt );
    
    if( PLOT )
        % [INFO] ...
        fprintf( 'Plotting bounds...' );
        
        % --- Plot bounds
        plotbnds( bdb7(:, :, i) );
        txt = ['Robust Tracking  Bounds for p' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot();
    end
end

% [INFO] ...
fprintf( ACK );

%% Step 8: Intersection of QFT Bounds and Compatibility

% [INFO] ...
fprintf( 'Step 8:' );
fprintf( '\tGrouping bounds...' );

% --- Grouping bounds
for i=1:width(P)
    bdb( :, :, i ) = grpbnds( bdb1(:,:,i), bdb2(:,:,i), bdb7(:,:,i) );
    if( PLOT )
        plotbnds( bdb( :, :, i ) );
        txt = ['All Bounds for p' num2str(i) num2str(i) '(s)' ];
        title( txt );
    end
end

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

for i=1:width(P)    
    % --- Find bound intersections
    ubdb( :, :, i ) = sectbnds( bdb( :, :, i ) );
    if( PLOT )
        plotbnds( ubdb( :, :, i ) );
        txt = ['Intersection of Bounds for p' num2str(i) num2str(i) '(s)' ];
        title( txt );
    end
end

% [INFO] ...
fprintf( ACK );

%% Step 9: Synthesize Feedback Controller G(s)

% [INFO] ...
fprintf( 'Step 9:' );
fprintf( '\tSynthesize G(s)...' );

% --- Directory where QFT generated controllers are stored
src = './controllerDesigns/';

for i=1:width(P)
    % --- Controller, G(s)
    G_file  = [ src 'g' num2str(i) num2str(i) '_i.shp' ];
    if( isfile(G_file) )
        g_ii( :, :, i ) = getqft( G_file );
    else
        if( i == 1 )
            num = -0.0006 .* [ 1/3e-5, 1];      % Numerator
            den = [ 1, 0 ];                     % Denominator
        elseif( i == 2 )
            num = -1.5 .* [ 1/4.5e-5, 1];       % Numerator
            den = [ 1, 0 ];                     % Denominator
        end
        
        % Construct controller TF
        g_ii( :, :, i ) = tf( num, den );
    end
end

% Define a frequency array for loop shaping
wl = logspace( log10(w(1)), log10(w(end)), 2048 );

% --- Loop over plants and design the controller
for i=1:width(P)
    L0(:, :, i) = P( i, i, nompt );
    L0(:, :, i).ioDelay = 0; % no delay
    lpshape( wl, ubdb(:, :, i), L0(:, :, i), g_ii( :, :, i ) );
    qpause;
end

% [INFO] ...
fprintf( ACK );


%% Step 10: Synthesize Prefitler F(s)

% [INFO] ...
fprintf( 'Step 10:' );
fprintf( '\tSynthesize F(s)...' );

for i=1:width(P)
    % --- Pre-filter, F(s)
    F_file  = [ src 'f' num2str(i) num2str(i) '_i.fsh' ];
    if( isfile(F_file) )
        f_ii( :, :, i ) = getqft( F_file );
    else
        if( i == 1 )
            num = 1;                            % Numerator
            den = [ 1/3.2e-5, 1 ];              % Denominator
        elseif( i == 2 )
            num = 1;                            % Numerator
            den = [ 1/3.2e-5, 1 ];              % Denominator
        end
        
        % Construct controller TF
        f_ii( :, :, i ) = tf( num, den );
    end
end

WW = logspace( log10( omega_6(1) ), ...         % Refine frequency array
               log10( omega_6(end) ), 1024 );

% --- Loop over plants and design the pre-filter
for i=1:width(P)
    PP = P( i, i, nompt );                      % Extract plant
    GG = g_ii( :, :, i );                       % Extract controller
    FF = f_ii( :, :, i );                       % Extract pre-filter

    % Loopshape
    pfshape( 7, WW, del_6, PP, [], GG, [], FF );
    qpause;
end

% [INFO] ...
fprintf( ACK );

%% Step 11-13: ANALYSIS

% [INFO] ...
fprintf( 'Steps 11-13:' );
fprintf( '\tRun Analysis...' );

for i=1:width(P)
    PP = P( i, i, nompt );                      % Extract plant
    GG = g_ii( :, :, i );                       % Extract controller
    FF = f_ii( :, :, i );                       % Extract pre-filter
    
    fprintf( "Stability Margins Specification\n" );
    fprintf( '\t> chksiso(1, wl, del_1, p_%i%i, [], g_%i%i, [], f_%i%i)\n', ...
              i, i, i, i, i, i)
    chksiso( 1, wl, del_1, PP, [], GG, [], FF );
    % [INFO] ...
    fprintf( "\t\t> " ); fprintf( ACK );
    
    fprintf( "Sensitivity Reduction Specification\n" );
    fprintf( '\t> chksiso(2, wl, del_3, p_%i%i, [], g_%i%i, [], f_%i%i)\n', ...
              i, i, i, i, i, i)
    ind = wl <= max(omega_3);
    chksiso( 2, wl(ind), del_3, PP, [], GG, [], FF );
    % [INFO] ...
    fprintf( "\t\t> " ); fprintf( ACK );

    fprintf( "Input Disturbance Rejection Specification\n" );
    fprintf( '\t> chksiso(7, wl, del_6, p_%i%i, [], g_%i%i, [], f_%i%i)\n', ...
              i, i, i, i, i, i)
    ind = find(wl <= max(omega_6));
    chksiso( 7, wl(ind), del_6, PP, [], GG, [], FF );
    % [INFO] ...
    fprintf( "\t\t> " ); fprintf( ACK );
end

% % [INFO] ...
% fprintf( ACK );

%% Check system/controller against Nyquist stability guidelines

% --- NOTE:
%   * Adding a zero corresponds to a +ve phase gain of +45deg / decade
%   * Adding a pole corresponds to a -ve phase drop of -45deg / decade
%
%   * Adding a  differentiator  shifts initial phase by +90deg
%   * Adding an integrator      shifts initial phase by -90deg
%
%   * For complex roots, phase gain/drop is +/-90deg

%{
% Open-loop TF
T_OL = P0*G;
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

%}

%% Check plant against Nyquist stability guidelines

%{
output = nyquistStability( P0 );

if( PLOT )
    figure();  rlocus( P0 ); grid on;
    figure(); nichols( P0 ); grid on;
    figure(); nyquist( P0 );
end
%}

%% MISC. TEMPORARY OPERATIONS

%{
clc;
% Open-loop TF
T_OL = P0*G;
% Closed-loop TF
T_CL = T_OL/(1+T_OL);
fprintf( "\n-> G(s)\n" ); nyquistStability( tf(G), false )
fprintf( "\n-> P(s)\n" ); nyquistStability( P0, false )
fprintf( "\n-> L(s)\n" ); nyquistStability( T_OL, false )
%}
