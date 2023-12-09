% Horizontal Wind Turbine QFT control design
%   Based on Dr. Mario Garcia-Sanz's book
%       "Robust Control Engineering: Practical QFT Solutions"
%       Example 8.1 - Pg. 187
%
%   Example used: 2x2 MIMO system
%       Method 2 MIMO Controller Design
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Sep. 30th, 2023
%
% CHANGELOG :
%   Sep. 30th, 2023
%       - Initial script
%

%% Setup environment
clear variables;
close all;
clc;
set( groot, 'defaultLineLineWidth', 1.5 );	% Set default line width of plots
set( 0, 'DefaultAxesFontSize', 12 );        % Set a default axes font size

% Change default interpreter (affects text, labels, etc...)
set( 0, 'DefaultTextInterpreter', 'latex' );
set( 0, 'DefaultLegendInterpreter', 'latex' );
set( 0, 'DefaultAxesTickLabelInterpreter', 'latex' );

format compact;
fontsize = 12;
tic

%% Flags/Constants

% --- Figure counter
CNTR = 1;                                   % Figure handle counter

% --- Enable/disable plotting figures
PLOT = true;                                %#ok<NASGU> If true, plot figures!
% PLOT = false;                               % COMMENT OUT TO PLOT FIGURES

% --- Enable/disable printing figures
PRNT = true;                                %#ok<NASGU>
PRNT = false;                               % COMMENT OUT TO PRINT FIGURES

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
if( isunix )
    ctrlSrc = [ filesep ctrlSrc ];
end

% --- Directory where QFT generated controllers are stored
dirG = fullfile( ctrlSrc, 'R3', 'G' );
dirF = fullfile( ctrlSrc, 'R3', 'F' );


%% Read A, B, C, D matrices from linearized model
data_dir    = './data/';
% name_mdl    = 'Simulink_Linearized_Model.mat';
% name_mdl    = 'Simulink_Linearized_Model_REV002.mat';
name_mdl    = 'Simulink_Linearized_Model_REV003.mat';
stateSpace  = load( [data_dir name_mdl ] );
% stateSpace  = stateSpace.Simulink_Linearization_Model_Timed_Based_Linearization;
stateSpace  = stateSpace.Simulink_Linearization_Model_FMIKit_Timed_Based_Linearization;

% % --- Get number of states
% [nStates, ~, ~] = size( stateSpace.A );
% --- Extract matrices from the one, big ABCD matrix
for i = 1:length( stateSpace )
    A_full(:,:,i) = stateSpace(i).a;
    B_full(:,:,i) = stateSpace(i).b;
    C_full(:,:,i) = stateSpace(i).c;
    D_full(:,:,i) = stateSpace(i).d;
end

% --- Choose the first matrix entry as the nominal
A = A_full( :, :, 1 );
B = B_full( :, :, 1 );
C = C_full( :, :, 1 );
D = D_full( :, :, 1 );

% --- Generate state-space model
% States and inputs names
% stateNames  = [ "phi"           , "omega"           , ...
%                 "blade120_phi"  , "blade120_omega"  , ...
%                 "blade0_phi"    , "blade0_omega"    , ...
%                 "blade240_phi"  , "blade240_omega"  , ...
%                 "V_{wind_x}"    , "V_{wind_y}"      , "V_{wind_z}"   ];
inputNames  = [ "u_{CPC}"       , ...
                "u_{IPC_{120}}" , "u_{IPC_{0}}"     , "u_{IPC_{240}}" ];
outputNames = [ "\omega_{rot}"  , ...
                "BRBM_{120}"    , "BRBM_{0}"        , "BRBM_{240}"];
% State-space model
sys         = ss( A, B, C, D                , ... $'StateName' , stateNames  , ...
                  'InputName' , inputNames  , ...
                  'OutputName', outputNames );
% --- Generate TF from SS model
TF = tf( sys );

% --- Generate TF using theory
syms s;
I = eye( size(A) );
P_manual = C*(s*I - A)^-1*B + D;


%% Actuator TF
wn_act = 0.1885;                % In rad/s
M_act = tf( 1, [1/wn_act 1] );  % Actuator dynamics TF


%% Step 1: Plant Modeling & Uncertainty

% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 300 ) == SISO with 300 TFs
%
% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 300 ) == SISO with 300 TFs
%

[~, ~, n_Plants] = size( A_full ) ;                 % Number of plants
[~, n_Inputs, ~] = size( B_full ) ;                 % Number of inputs
[n_Outputs, ~,~] = size( D_full ) ;                 % Number of outputs

P   = tf( zeros(n_Outputs,n_Inputs,n_Plants) );     % Pre-allocate memory

Pdiag   = tf( zeros(n_Outputs,n_Inputs,n_Plants) ); % Pre-allocate memory
Pinv    = tf( zeros(n_Outputs,n_Inputs,n_Plants) ); % Pre-allocate memory
PinvPdiag = tf(zeros(n_Outputs,n_Inputs,n_Plants)); % Pre-allocate memory

gain_PinvPdiag = zeros( size(PinvPdiag) );          % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

NDX = 1;                                            % Plant counter
for variation = 1:n_Plants                          % Loop over variations
    
    % --- Generate grided TF from grided SS model
    A_g = A_full(:, :, variation);
    B_g = B_full(:, :, variation);
    C_g = C_full(:, :, variation);
    D_g = D_full(:, :, variation);

    sys_g = ss( A_g, B_g, C_g, D_g ); sys_g_original = sys_g;
    [sys_g, prescale_info] = prescale( sys_g );
    TF_g = tf( sys_g );

    % --- Here we create the plant TF
    % --- 1st row
    p11_NDX        = TF_g( 1, 1 );
    p12_NDX        = TF_g( 1, 2 );
    p13_NDX        = TF_g( 1, 3 );
    p14_NDX        = TF_g( 1, 4 );
    
    % --- 2nd row
    p21_NDX        = TF_g( 2, 1 );
    p22_NDX        = TF_g( 2, 2 );
    p23_NDX        = TF_g( 2, 3 );
    p24_NDX        = TF_g( 2, 4 );

    % --- 3rd row
    p31_NDX        = TF_g( 3, 1 );
    p32_NDX        = TF_g( 3, 2 );
    p33_NDX        = TF_g( 3, 3 );
    p34_NDX        = TF_g( 3, 4 );

    % --- 4th row
    p41_NDX        = TF_g( 4, 1 );
    p42_NDX        = TF_g( 4, 2 );
    p43_NDX        = TF_g( 4, 3 );
    p44_NDX        = TF_g( 4, 4 );
    
    % --- Place them all in one big matrix
    % ***NOTE:
    %   P( 1, 1,  1, : ) ==> p11 / 1st  variation (i.e. p11(:,:,1))
    %   P( 2, 1, 65, : ) ==> p21 / 65th variation (i.e. p21(:,:,65))
    P(:, :, NDX) = [ p11_NDX, p12_NDX, p13_NDX, p14_NDX ;
                     p21_NDX, p22_NDX, p23_NDX, p24_NDX ;
                     p31_NDX, p32_NDX, p33_NDX, p34_NDX ;
                     p41_NDX, p42_NDX, p43_NDX, p44_NDX ];

%     % --- Generate modified plant by inclusion of actuator dynamics
%     P(:, :, NDX) = P(:, :, NDX) .* M_act;

    % --- Generate diagonal matrix, Pdiag(s)
    Pdiag(:, :, NDX) = [ P(1,1,NDX),    0       ,   0        ,  0       ;
                            0      , P(2,2,NDX) ,   0        ,  0       ; 
                            0      ,    0       , P(3,3,NDX) ,  0       ;
                            0      ,    0       ,   0        , P(4,4,NDX) ];

    % --- Generate inverted matrix, Pinv(s)
    Pinv(:, :, NDX) = minreal( inv(P(:, :, NDX)) );

    % --- Generate temporary G_α(s) = Pinv(s) * Pdiag(S)
    PinvPdiag(:, :, NDX) = Pinv(:, :, NDX) * Pdiag(:, :, NDX);
    PinvPdiag(:, :, NDX) = minreal( PinvPdiag(:, :, NDX), 0.01 );

    % --- Get the DC gain
    gain_PinvPdiag(:, :, NDX) = dcgain( PinvPdiag(:, :, NDX) );

    NDX = NDX + 1;                      % Increment counter
    if( ~rem(NDX, 50) ); disp(NDX); end
end


% % --- EXTRA STEP: modify Pw(s) as per the problem requirement
% % Add low-ass filter
% f_LP      = tf( 1, [1/(1.1e-4)^2, 2/(1.1e-4), 1 ] );
% % Generate modified plants
% P = Pw.*f_LP;

% [INFO] ...
fprintf( ACK );

% Cleanup
clearvars min_* max_* grid_*


%% Step 2: The Nominal Plant

% [INFO] ...
fprintf( 'Step 2:' );
fprintf( '\tComputing nominal plant...' );

% --- Generate nominal plant TF
%   Any one of the models above can be used as the nominal plant.
%   We just happened to chose this one.
%
% P0 = TF;                        % Nominal Transfer Function
% Get average
A0 = mean( A_full, 3 );     B0 = mean( B_full, 3 );
C0 = mean( C_full, 3 );     D0 = mean( D_full, 3 );
P0 = prescale( ss(A0, B0, C0, D0) );
P0 = tf( P0 );              % Nominal Transfer Function
P(:, :, end+1) = P0;
nompt = length(P);

% --- Define nominal plant case (recall, P(:,:,1) corresponds to P0)
% nompt = 8;

% --- Get total plants size
% [x-dim, y-dim, z-dim] = [nrowsP, ncolsP, nvarsP]
[nrowsP, ncolsP, nvarsP] = size( P );

% --- Cleanup plants transfer function by removing values below 1e-08 and
% minreal of 0.01
% P = numerical_cleanup( P, 1e-08, 0.01 );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
ww = logspace( log10(0.10), log10(50.0), 2048 );
[p0, theta0] = bode( P0, ww );
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    bode( P0, '-', P0, '.r', ww(1:32:end) ); grid on;
    make_nice_plot( PRNT, './figs', 'bode_plot' );
end

% Cleanup
clearvars A* B* C* D* var* *_NDX -except A B C D ACK CNTR


%% Step 3: QFT Template

% [INFO] ...
fprintf( 'Step 3:' );
fprintf( '\tPlotting QFT templates...' );

% --- Working frequencies
% w = [ 0.5 1.0 5.0 6.0 7.0 8.0 9.0 10.0 ];
w = [ 0.10 0.25 0.50 0.75 1.00 2.50 ...
      5.00 6.00 7.00 8.00 9.00 10.0 ...
      15.0 20.0 25.0 30.0 35.0 40.0 50.0 ];

if( PLOT )
    % --- Plot QFT templates
    for ROW = 1:nrowsP
        for COL = 1:ncolsP
            plottmpl( w, P(ROW, COL, :, :), nompt );
    
            % --- Change legend position
            hLegend = findobj( gcf, 'Type', 'Legend' ); % Get legend property
            set( hLegend, 'location', 'southeast' );    % Access and change location
            
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
P0_0    = dcgain( P0 );
[U,S,V] = svd( P0_0 );
P0_0inv = V/S*U.';
Lambda_0 = P0_0 .* P0_0inv.';

% === NEED REVISION
% === FOR NOW, INSPECT BY HAND!
% --- Store maximum and minimum gains allowable for the MIMO system
% maxGain = -inf; dirMaxGain = 0;
% minGain = +inf; dirMinGain = 0;
% for ii = 1:length(S)
%     tempMax = max( S(:,ii), [], "all" );
%     tempMin = min( S(:,ii), [], "all" );
% 
%     if( tempMax > maxGain )
%         maxGain = tempMax;
%         dirMaxGain = ii;
%     end
% 
%     if( tempMin < minGain )
%         minGain = tempMin;
%         dirMinGain = ii;
%     end
% end

% --- RGA for s=jw=inf
P0_inf = freqresp( P0, 1e16 );
P0_inf = abs( P0_inf );                 % Make sure we get sensical numbers
[U,S,V] = svd( P0_inf );
P0_infinv = V/S*U.';
Lambda_inf = P0_inf .* P0_infinv.';

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
    
    % --- Cast numbers less than zero (~=1e-8 being close enough) as -inf
    for ii = 1:length( Lambda_COL )
        if( Lambda_COL(ii) < 1e-8 )
            Lambda_COL(ii) = -inf;
        end
    end
    
    % --- Carry on finding closest number to one that is non-negative
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
omega_1 = [ 0.10 0.25 0.50 0.75 1.00 2.50 ...
            5.00 6.00 7.00 8.00 9.00 10.0 ];

% Restriction (for p_ii, i=1,2)
% W_s         = 1.66;
% W_s         = 1.46;
W_s         = 1.08;
% W_s         = 1.01;
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
% omega_3 = [ 5e-2 7.5e-2 1e-1 ];
omega_3 = [ 0.10 0.25 0.50 0.75 1.00 ];

% Restriction
a_d     = 2;
num     = [ 1/a_d   , 0 ];
den     = [ 1/a_d   , 1 ];
del_3   = tf( num, den );

% --- Plot bounds
if( PLOT )
    figure( CNTR ); CNTR = CNTR + 1;
    w_del_3 = logspace( log10(w(1)), log10(w(end)));
    [mag, ~] = bode( del_3, w_del_3 );
    mag_dB = db( squeeze(mag) );

    semilogx( w_del_3, mag_dB ); grid on;
    
    txt = ["Sensitivity Specification"];
    title( txt );
    make_nice_plot( PRNT, './figs', txt );
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
% omega_6 = [ 5e-2 7.5e-2 1e-1 2.5e-1 ];
omega_6 = [ 0.10 0.25 0.50 0.75 1.00 2.50 5.00 ];


% Restriction
% -----------
% Upper bound
% -----------
a_U = 1.8; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.00;
num = [ conv([1/a_U 1], [0 1+eps_U]) ];
den = [ (1/wn)^2 (2*zeta/wn) 1 ];
del_6_U = tf( num, den );
% -----------
% Lower bound
% -----------
a_L = 4.4; eps_L = 0.00;
num = 1-eps_L;
den = [ conv([1/a_L 1], [1/a_L 1]) ];
del_6_L = tf( num, den );
% -----------
% Tracking weight
% -----------
del_6 = [ del_6_U  ;
          del_6_L ];

if( PLOT )
    % --- PLOT step response of del_6_U(s) and del_6_L(s)
    figure( CNTR ); CNTR = CNTR + 1;
    stepplot( del_6_U, del_6_L ); grid on;
    title( "Reference Tracking Specification" );
    make_nice_plot();
end

% [INFO] ...
fprintf( ACK );

%% Step 9.1: Synthesize Feedback Controller G_α(s)

% --- The fully populated matrix controller G(s) is composed of two
% matrices: G(s) = G_alpha(s)*G_beta(s)
%
%   Let α = alpha and β = beta. Then,
%                      _                     _     _                     _
%                     |  g_11α  g_12α  g_13α  |   |  g_11β    0       0   |
%       G = G_α*G_β = |  g_21α  g_22α  g_23α  | * |    0    g_22β     0   |
%                     |_ g_31α  g_32α  g_33α _|   |_   0      0    g_33β _|
%
%   The main objective of the pre-compensator is to diagnolize the plant
%   P(s) as much as possible. Therefore, the expression used to calculate
%   G_a(s) is,
%              _                     _ 
%             |  g_11α  g_12α  g_13α  |
%       G_α = |  g_21α  g_22α  g_23α  | = P(s)^-1*P_diag(s)     (Eqn. 1)
%             |_ g_31α  g_32α  g_33α _|
%
%              _                     _     _                  _
%             |  p'_11  p'_21  p'_31  |   |  p_11   0     0    |
%           = |  p'_21  p'_22  p'_23  | * |   0    p_21   0    |
%             |_ p'_31  p'_32  p'_33 _|   |_  0     0    p_33 _|
%
%   Where p'_ij corresponds to ij-th element of the inverted P(s) matrix.
%
%
%   The role of the controller G_α is not to achieve an exact decoupling,
% but to ease the design of G_β. That is, to reduce the amount of feedback
% needed to achieve the robust performance specifications. Besides, this
% approach allows modifying, when necessary, the final form of the
% controller G_α in Eqn. 1 so that:
%
%   -   No RHP or imaginary axis pole-zero cancellation occurs between P
%       and G_α or its elements,
%   -   No RHP transmission elements (Smith-McMillan) are introduced by
%       the controller G_α,
%   -   The relative difference of the number of poles and zeros in each
%       element of the G_α matrix is the same as in Eqn. 1 in order to
%       ease the design of the G_β controller,
%   -   The RGA of the system is improved, looking for positive and
%       close-to-one diagonal elements λ_ii of the RGA matrix. That is,
%       the pre-compensator G_α decouples the system to some extent, which
%       is its main goal.
%

% [INFO] ...
fprintf( 'Step 9:' );
fprintf( '\tSynthesize G(s)\n' );

% --- Working frequency range
wl = logspace( log10(w(1)), log10(w(end)), 2048 );

% --------------------------------------------------
% ----              Generate G_α(s)             ----
% --------------------------------------------------

% [INFO] ...
fprintf( '\t> Computing G_alpha(s)...' );

% % --- Generate diagonal matrix
% P_diag  = tf( zeros(size(P)) );             % Pre-allocate memory
% for ii  = 1:width( P )
%     P_diag( ii, ii, :, : )  = P( ii, ii, :, : );
% end
% 
% % --- Compute the gain of all elements in P^1(s) * P_diag(s)
% TOL = 0.01;                                 % Tolerance for cancellation
% Pinv      = inv( P );                       % Compute P^-1
% PinvPdiag = minreal( Pinv * P_diag, TOL );  % Compute P^-1*P_diag
% 
% gain_PinvPdiag = zeros( size(P) );          % Pre-allocate memory
% for ROW = 1:width( PinvPdiag )              % Loop over ROWS
%     for COL = 1:width( PinvPdiag )          % Loop over COLS
%         for NDX = 1:n_Plants                % Loop over variations
% 
%             % Get the n-th plant
%             nth_Plant = PinvPdiag(ROW, COL, NDX, :);
%             % Compute DC gain
%             kP = dcgain( nth_Plant );
%             % Store in a matrix
%             gain_PinvPdiag(ROW, COL, NDX, :) = kP;
% 
%         end  
%     end
% end

% --- Compute the mean value
NROWS = width( gain_PinvPdiag );
NCOLS = width( gain_PinvPdiag );
meanGain_PinvPdiag = zeros( NROWS, NCOLS ); % Pre-allocate memory
for ROW = 1:NROWS                           % Loop over ROWS
    for COL = 1:NCOLS                       % Loop over COLS

        meanGain_PinvPdiag(ROW, COL) = mean( gain_PinvPdiag(ROW, COL, :) );

    end
end

% --- Lastly, construct initial G_α(s) controller based on the
%   mean value obtained.
%
%   ***NOTE:
%       This is NOT necessarily the final form of G_α(s), as we may
%   need/want to tweak it to avoid certain frequencies for instance. 
%

err         = inf.* ones( size(P0) );
newErr      = zeros( size(P0) );
nPinvPdiag  = zeros( size(P0) );
G_alpha     = tf( zeros(size(nPinvPdiag)) );

for ROW = 1:width( PinvPdiag )              % Loop over ROWS
    for COL = 1:width( PinvPdiag )          % Loop over COLS
        for NDX = 1:n_Plants                % Loop over variations
            newErr(ROW, COL) = abs( meanGain_PinvPdiag(ROW, COL) - gain_PinvPdiag( ROW, COL, NDX) );
            if( newErr(ROW, COL) <= err(ROW, COL) )
                nPinvPdiag(ROW, COL) = NDX;
                err(ROW, COL) = newErr(ROW, COL);
            end
        end
    end
end

% --- g_alpha_ij(s) is defined here
for ROW = 1:width( PinvPdiag )              % Loop over ROWS
    for COL = 1:width( PinvPdiag )          % Loop over COLS

        g_ij = PinvPdiag( ROW, COL, nPinvPdiag(ROW, COL) );
        G_alpha( ROW, COL ) = minreal( g_ij, 0.1 );

        % figure(); bode(g_ij); hold on; grid on;
        % bode(minreal( g_ij, 0.1 )); make_nice_plot();

    end
end

TOL = 0.1;
% TOL = 1.0;
% -----------
% --- 1ST ROW
% -----------
% g11_a = minreal( G_alpha(1, 1), TOL );                      % Extract controller
% controlSystemDesigner( 'bode', 1, g11_a );                  % Loop-shape
% qpause;
g11_a = tf( [-1.5860e+03 -1.3947e+06 -2.9549e+08 -2.9872e+10 -1.2753e+12 -1.4708e+13 -8.6190e+11], ...
            [1 697.3444 1.6325e+05 8.1728e+05 1.4976e+05 1.3491e+04 366.4242] );     % Updated Tuned controller


% g12_a = minreal( G_alpha(1, 2), TOL );                      % Extract controller
% controlSystemDesigner( 'bode', 1, g12_a );                  % Loop-shape
% qpause;
g12_a = tf( [-8.3255e+08 -1.2702e+11 -5.6845e+12 -6.4799e+13 -2.9971e+12 -1.7625e+08], ...
            [1 692.4138 1.5984e+05 3.1283e+04 120.4510 0.2304] );    % Updated Tuned controller


g13_a = minreal( G_alpha(1, 3), TOL );                      % Extract controller
% controlSystemDesigner( 'bode', 1, g13_a );                  % Loop-shape
% qpause;
% g13_a = tf( [2.08e06 4.571e08], [1 3435 -343.2 858.7] );    % Updated Tuned controller


g14_a = minreal( G_alpha(1, 4), TOL );                      % Extract controller
% controlSystemDesigner( 'bode', 1, g14_a );                  % Loop-shape
% qpause;
% g14_a = tf( [2.428e05 1.214e08], [1 1500 -224.4 843.8] );   % Updated Tuned controller


% -----------
% --- 2ND ROW
% -----------
g21_a = minreal( G_alpha(2, 1), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g21_a );  % Loop-shape
% qpause;
% g21_a = tf( [-600 -1200 -600], [1 8.75 6] );  % Updated Tuned Controller
% g21_a = tf( 100.8, [1 10] );                  % Updated Tuned controller


g22_a = minreal( G_alpha(2, 2), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g22_a );  % Loop-shape
% qpause;
% g22_a = tf( [187.5 37.5], [1 0.5 -1.5] );     % Updated Tuned Controller
% g22_a = tf( [3.6414e+03 9.4036e+06], [1 -1.5041 314.2805] );      % Updated Tuned controller***CHECK


g23_a = minreal( G_alpha(2, 3), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g23_a );  % Loop-shape
% qpause;
% g23_a = tf( [3.6414e+03 9.4036e+06], [1 -1.5041 314.2805] );      % Updated Tuned controller***CHECK


g24_a = minreal( G_alpha(2, 4), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g24_a );  % Loop-shape
% qpause;
% g24_a = tf( [4.3875e+03 1.6155e+07 3.5758e+09], [1 -159.7504 -223.9247 -8.8853e+04] );      % Updated Tuned controller***CHECK


% -----------
% --- 3RD ROW
% -----------
g31_a = minreal( G_alpha(3, 1), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g31_a );  % Loop-shape
% qpause;
% g31_a = tf( 39.85, [1 3.953] );               % Updated Tuned controller


g32_a = minreal( G_alpha(3, 2), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g32_a );  % Loop-shape
% qpause;
% g32_a = tf( [3.6414e+03 9.4036e+06], [1 -1.5041 314.2805] );      % Updated Tuned controller***CHECK


g33_a = minreal( G_alpha(3, 3), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g33_a );  % Loop-shape
% qpause;
% g33_a = tf( 0.2216, [1 0.00015] );          % Updated controller
% g33_a = tf( 100.8, [1 10] );                  % Updated Tuned controller


g34_a = minreal( G_alpha(3, 4), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g34_a );  % Loop-shape
% qpause;
% g34_a = tf( [3.6414e+03 9.4036e+06], [1 -1.5041 314.2805] );      % Updated Tuned controller***CHECK


% -----------
% --- 4TH ROW
% -----------
g41_a = minreal( G_alpha(4, 1), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g41_a );  % Loop-shape
% qpause;
% g41_a = tf( 39.85, [1 3.953] );               % Updated Tuned controller


g42_a = minreal( G_alpha(4, 2), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g42_a );  % Loop-shape
% qpause;
% g42_a = tf( 100.8, [1 10] );                  % Updated Tuned controller


g43_a = minreal( G_alpha(4, 3), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g43_a );  % Loop-shape
% qpause;
% g43_a = tf( 100.8, [1 10] );                  % Updated Tuned controller


g44_a = minreal( G_alpha(4, 4), TOL );      % Extract controller
% controlSystemDesigner( 'bode', 1, g44_a );  % Loop-shape
% qpause;
% g44_a = tf( [3.6414e+03 9.4036e+06], [1 -1.5041 314.2805] );      % Updated Tuned controller***CHECK


G_alpha = [ g11_a, g12_a, g13_a, g14_a ;
            g21_a, g22_a, g23_a, g24_a ;
            g31_a, g32_a, g33_a, g34_a ;
            g41_a, g32_a, g33_a, g44_a ];

% [INFO] ...
fprintf( ACK );

% --- Plot to visualize
if( ~PLOT )
    for ROW = 1:width( gain_PinvPdiag )         % Loop over ROWS
        for COL = 1:width( gain_PinvPdiag )     % Loop over COLS
            figure(); bode( PinvPdiag(ROW, COL, 1:1:end, :), wl );  grid on;
            hold on ; bode( G_alpha( ROW, COL ), wl(1:16:end), 'r*' );
            bode( G_alpha( ROW, COL ), wl(1:16:end), 'r--' );
            
            text_1 = [ 'p*_{' num2str(ROW) num2str(COL) '}(s)' ];
            text_2 = [  'p_{' num2str(COL) num2str(COL) '}(s)' ];
            text_3 = [  'g_{' num2str(ROW) num2str(COL) '}(s)' ];
            title( [text_1 ' \times ' text_2 ' and ' text_3] );

            txt = ['pstar_' num2str(ROW) num2str(COL) ' x ' ...
                   'p_' num2str(COL) num2str(COL) ' and ' ...
                   'g_' num2str(ROW) num2str(COL) ];
            make_nice_plot( PRNT, './figs', txt );
        end
    end
end


% --- As we can see, the initial G_α(s) controller based on the
%   mean value works great for low frequencies. However, we need it to
%   then filter out the dynamics before the nmp zero at −2 × 10–4
%   rad/s
%
%   Let's use the Control System Designer Toolbox to loopshape the
%   G_α(s) = g_α_ij controller and reiterate over this process
%

% Cleanup
% clearvars g* -except gain_PinvPdiag

%% ================================================
%  ===== STOPPED HERE. CONTINUE EDITING BELOW =====
%  ================================================

%% Step 9.1: Synthesize Feedback Controller G_β(s)

% --------------------------------------------------
% ----              Generate G_β(s)             ----
% --------------------------------------------------
%
%   Recall, G_β(s) is the diagonal controller matrix.
%              _                     _
%             |  g_11β    0       0   |
%       G_β = |    0    g_22β     0   |
%             |_   0      0    g_33β _|

% [INFO] ...
fprintf( '\t> Computing G_beta(s)...' );

% --- Start by computing the extended matrix Px = P*G_α
%
%   Recall,
%
%       In addition, the plant matrix P(s), the corresponding inverse
%   P(s)^–1, and the diagonal P_diag(s) are selected so that the expression
%   of the extended matrix Px = P*G_α presents the closest form to a
%   diagonal matrix, nulling the off-diagonal terms as much as possible.
%

TOL     = 0.1;
Px      = minreal( P * G_alpha, TOL );              % Extended matrix
Pxinv   = minreal( inv(Px), TOL );                  % Invert extended matrix
nulledP = abs( freqresp(Px(:,:,end), 0.5, 'rad/s') );

% Pxsvd = tf( zeros(size(Px)) );
% tic;
% for ii = 1:length(Px)
%     for ROW = 1:width( Pxsvd )                      % Loop over ROWS
%         for COL = 1:width( Pxsvd )                  % Loop over COLS
%             [U,S,V] = svd( Px(ROW,COL,ii,:) );
%             Pxsvd(ROW,COL,ii,:) = minreal( V/S*U', TOL );
%         end
%     end
% end
% toc;

% --- Sequential desgin (loop-shape) gbeta_ii(s), where
%
%       > g_ij(s)  = 0 for i != j
%       > g_ij(s) != 0 for i  = j
%
%   Furthermore, recall that the loop L_ii(s) is defined as:
%
%       > L_ii(s) = qx_ii(s) * gbeta_ii(s)
%
%   Where qx_ii(s) = 1/px_ii(s) = big expression
%

% --- Loopshape g11_b(s) controller over L11(s) = qx11(s) * g11_b(s)
%

% qx11 = tf( zeros(NROWS, NCOLS) );                   % Pre-allocate memory
qx11( 1, 1, : ) = 1/Pxinv( 1, 1, : );             % Initialize



% --- Directory where QFT generated controllers are stored
src = './controllerDesigns/';

% --- Controller, G(s)
% --- Controller, G(s)
G_name  = ['g' num2str(i) num2str(i) '.shp'];
G_file  = fullfile( dirG, G_name );
if( isfile(G_file) )
    fprintf( "\tController %s found. Loading from file.\n", G_name );
    g_ii( :, :, i ) = getqft( G_file );
else
    syms s;
    num = (0.9) .* sym2poly(   (s/0.24 + 1) );      % Numerator
    den =          sym2poly( s*(s/40.0 + 1) );      % Denominator
    clear s;
    
    % Construct controller TF
    g11_b = tf( num, den );                         % Eq.(CS3.25)
end

Q = { qx11 };

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
for i=1:1
    p_ii = Q{ i };
    bdb1(:, :, i) = sisobnds( spec, omega_1, del_1, p_ii, [], nompt );
    % R = 0; bdb1 = sisobnds( spec, omega_1, del_1, P, R, nompt );

    if( PLOT )
        % [INFO] ...
        fprintf( 'Plotting bounds...' );

        % --- Plot bounds
        plotbnds( bdb1(:, :, i) );
        txt = ['Robust Stability Bounds for qx' num2str(i) num2str(i) '(s)' ];
        title( txt );
        % xlim( [-360 0] ); ylim( [-10 30] );
        make_nice_plot( PRNT, './figs', txt );
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
for i=1:1
    p_ii = Q{ i };
    bdb2(:, :, i) = sisobnds( spec, omega_3, del_3, p_ii, [], nompt );

    if( PLOT )
        % [INFO] ...
        fprintf( 'Plotting bounds...' );

        % --- Plot bounds
        plotbnds( bdb2(:, :, i) );
        txt = ['Sensitivity Reduction Bounds for qx' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot( PRNT, './figs', txt );
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
for i=1:1
    p_ii = Q{ i };
    bdb7(:, :, i) = sisobnds( spec, omega_6, del_6, p_ii, [], nompt );

    if( PLOT )
        % [INFO] ...
        fprintf( 'Plotting bounds...' );

        % --- Plot bounds
        plotbnds( bdb7(:, :, i) );
        txt = ['Robust Tracking  Bounds for qx' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot( PRNT, './figs', txt );
    end
end

% [INFO] ...
fprintf( ACK );

%% Step 8: Intersection of QFT Bounds and Compatibility

% [INFO] ...
fprintf( 'Step 8:' );
fprintf( '\tGrouping bounds...' );

% --- Grouping bounds
for i=1:1
    bdb( :, :, i ) = grpbnds( bdb1(:,:,i), bdb2(:,:,i), bdb7(:,:,i) );
    if( PLOT )
        plotbnds( bdb( :, :, i ) );
        txt = ['All Bounds for qx' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot( PRNT, './figs', txt );
    end
end

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

for i=1:1 
    % --- Find bound intersections
    ubdb( :, :, i ) = sectbnds( bdb( :, :, i ) );
    if( PLOT )
        plotbnds( ubdb( :, :, i ) );
        txt = ['Intersection of Bounds for qx' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot( PRNT, './figs', txt );
    end
end

% [INFO] ...
fprintf( ACK );

%% Start loopshaping GUI
L11 = qx11( 1, 1, nompt );                          % Desired loop
L11.ioDelay = 0;                                    % No delay
lpshape( wl, ubdb(:, :, 1), L11, g11_b );
% qpause;

% --- Loopshape g22_b(s) controller over L22(s) = qx22(s) * g22_b(s)
px22 = tf( zeros(1, 1, length(Pxinv)) );
qx22 = tf( zeros(1, 1, length(Pxinv)) );
for ii = 2:2
    % --- Use sequential method
    gg = Pxinv(ii, ii, :) - ...
         ((Pxinv(ii, ii-1, :) * Pxinv(ii-1, ii, :)) / ...
         (Pxinv(ii-1, ii-1, :) + g11_b));
    px22( 1, 1, : ) = minreal( gg, 0.1 );
end
qx22( 1, 1, : ) = 1/px22( 1, 1, : );


% % --- Cleanup transfer function by removing values below 1e-10
% for ii = 1:length( px22 )
%     qx22( 1, 1, ii ) = 1/px22( 1, 1, ii );
% 
%     [n, d] = tfdata(px22( 1, 1, ii ));
%     n = cellfun(@(x) {x.*(x>1e-16)}, n);
%     d = cellfun(@(x) {x.*(x>1e-16)}, d);
%     qx22( 1, 1, ii ) = 1/tf(n, d);
% end

% --- Controller, G(s)
G_file  = [ src 'HX_g22_b.shp' ];
if( isfile(G_file) )
    g22_b = getqft( G_file );
else
    syms s;
    num = (0.4) .* sym2poly(   (s/0.10 + 1) );      % Numerator
    den =          sym2poly( s*(s/40.0 + 1) );      % Denominator
    clear s;
    
    % Construct controller TF
    g22_b = tf( num, den );                         % Eq.(CS3.29)
end

Q = {qx11; qx22};


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
for i=2:2
    p_ii = Q{ i };
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
for i=2:2
    p_ii = Q{ i };
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
for i=2:2
    p_ii = Q{ i };
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
for i=2:2
    bdb( :, :, i ) = grpbnds( bdb1(:,:,i), bdb2(:,:,i), bdb7(:,:,i) );
    if( PLOT )
        plotbnds( bdb( :, :, i ) );
        txt = ['All Bounds for p' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot();
    end
end

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

for i=2:2 
    % --- Find bound intersections
    ubdb( :, :, i ) = sectbnds( bdb( :, :, i ) );
    if( PLOT )
        plotbnds( ubdb( :, :, i ) );
        txt = ['Intersection of Bounds for p' num2str(i) num2str(i) '(s)' ];
        title( txt );
        make_nice_plot();
    end
end

% [INFO] ...
fprintf( ACK );

%% <<<<<<<<<<<<<<<<<<<<<<                >>>>>>>>>>>>>>>>>>>>>>

% Start loopshaping GUI
L22 = qx22( 1, 1, nompt );                          % Desired loop
L22.ioDelay = 0;                                    % No delay
lpshape( wl, ubdb(:, :, 2), L22, g22_b );
% qpause;


% [INFO] ...
fprintf( ACK );
toc

%% Step 10: Synthesize Prefitler F(s)

% --- Loopshape pre-filters
f11 = tf( [1/20 1], [1/2 1] ); % Eq.(8.184)
f22 = tf( [1/20 1], [1/2 1] ); % Eq.(8.188)


%% Miscellaneous tests


aaa = [ -0.01   , 10  , 5     ;
         5     , -0.01, 10    ;
         10     , 5   , -0.01 ];

for COL = 1:width(aaa)
    aaa_COL = aaa(:, COL);                      % Extract column
    
    % --- Cast numbers less than zero as -inf
    for ii = 1:length( aaa_COL )
        if( aaa_COL(ii) < 0 )
            aaa_COL(ii) = -inf;
        end
    end
    
    % --- Carry on finding closest number to one that is non-negative
    abs_val = abs(aaa_COL-VAL);
    [minValue, NDX] = min( abs_val );
    closestValue = aaa_COL( NDX );

    fprintf( "\t> ( u%i, y%i )\n", COL, NDX );
end

