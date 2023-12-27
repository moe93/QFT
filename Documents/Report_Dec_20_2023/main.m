% Comparison between CPC vs CPC+IPC performance for fatigue load reduction
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Nov. 11th, 2023
%
% CHANGELOG :
%   Dec. 20th, 2023
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


%% Flags/Constants

% --- Figure counter
CNTR = 1;                                   % Figure handle counter

% --- Enable/disable plotting figures
PLOT = true;                                %#ok<NASGU> If true, plot figures!
% PLOT = false;                               % COMMENT OUT TO PLOT FIGURES

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
scriptDir   = mfilename( 'fullpath' );
scriptParts = strsplit( scriptDir, filesep );
% Go up one level and generate new path
src = fullfile( scriptParts{1:end-1} );

% If on a UNIX machine (i.e. macOS, Ubuntu, etc...), fix path since
% strsplit() removes the leading '/' from the path.
if( isunix )
    src = [ filesep src ];
end

% --- Data and figures directories
data_dir = fullfile( src, 'data' );
figs_dir = fullfile( src, 'figures' );

% --- Clear workspace
clear script* src

%% Load data

% --- CPC
CPC = load( fullfile( data_dir, 'QFT_AzimuthVar_4x4IndepSISO_NoIPC.mat') );
% Unpack
t_CPC        	= CPC.data(:, 1);           % Simulation time       {  s  }
RotSpd_CPC   	= CPC.data(:, 2);           % Rotor angular speed   {rad/s}
BRBM_CPC        = CPC.data(:, 3);           % BRBM                  { N.m }
betaCPC_CPC     = rad2deg(CPC.data(:, 4));  % Collective pitch      { deg }
betaIPC_CPC     = rad2deg(CPC.data(:, 5));  % Individual pitch      { deg }
GenTrq_CPC      = -1.*CPC.data(:, 6);       % Generator torque      { N.m }
Power_CPC       = GenTrq_CPC.*RotSpd_CPC;   % Power                 {  W  }

% Store into a cell (i.e a matrix of vectors) for easier access
CPC_Data	    = { t_CPC, RotSpd_CPC, betaCPC_CPC, BRBM_CPC, ...
                    betaIPC_CPC, GenTrq_CPC, Power_CPC };

% --- CPC + IPC
IPC = load( fullfile( data_dir, 'QFT_AzimuthVar_4x4IndepSISO_Simulink_LOWER_GAIN.mat') );
% Unpack
t_IPC        	= IPC.data(:, 1);           % Simulation time       {  s  }
RotSpd_IPC   	= IPC.data(:, 2);           % Rotor angular speed   {rad/s}
BRBM_IPC        = IPC.data(:, 3);           % BRBM                  { N.m }
betaCPC_IPC     = rad2deg(IPC.data(:, 4));  % Collective pitch      { deg }
betaIPC_IPC     = rad2deg(IPC.data(:, 5));  % Individual pitch      { deg }
GenTrq_IPC      = +1.*IPC.data(:, 6);       % Generator torque      { N.m }
Power_IPC       = GenTrq_IPC.*RotSpd_IPC;   % Power                 {  W  }

% Store into a cell (i.e a matrix of vectors) for easier access
IPC_Data	    = { t_IPC, RotSpd_IPC, betaCPC_IPC, BRBM_IPC, ...
                    betaIPC_IPC, GenTrq_IPC, Power_IPC };

% --- Clear workspace
clear CPC *_CPC IPC *_IPC data_dir


%% Process data

% --- FFT BRBM
% We wish to get the FFT starting at time t=400s forwards.
% Since our sampling time is 0.01, we can get the index as follows:
NDX = 400/0.01;
[f_CPC, P_CPC] = do_fft(CPC_Data{4}(NDX:end), CPC_Data{1}(NDX:end), 0.01);
[f_IPC, P_IPC] = do_fft(IPC_Data{4}(NDX:end), IPC_Data{1}(NDX:end), 0.01);

figure( CNTR ); CNTR = CNTR + 1;
plot( f_CPC, P_CPC ); hold on;
plot( f_IPC, P_IPC ); grid on; hold off;
xlabel( 'Frequency, $\left( Hz \right)$' )
ylabel( 'Amplitude, $\left( N.m \right)$' )

xlim( [0.00 0.25] );
ylim( [0.00 6e+6] );

title( 'PSD of BRBM' );

% Add legend
legend( 'CPC', 'CPC+IPC', ...
        'NumColumns'    , 1             , ...
        'Location'      , 'northeast' );
make_nice_plot( PRNT, figs_dir, 'Updated_PSD' );

% Get percent difference
max_fft_CPC = max( P_CPC(2:end) );
max_fft_IPC = max( P_IPC(2:end) );
prct_diff_fft = abs(max_fft_CPC - max_fft_IPC)/max_fft_CPC

% --- PSD BRBM
% We wish to get the PSD starting at time t=400s forwards.
% Since our sampling time is 0.01, we can get the index as follows:
NDX = 400/0.01;
[f_CPC, P_CPC] = do_psd(CPC_Data{4}(NDX:end), CPC_Data{1}(NDX:end), 0.01);
[f_IPC, P_IPC] = do_psd(IPC_Data{4}(NDX:end), IPC_Data{1}(NDX:end), 0.01);

figure( CNTR ); CNTR = CNTR + 1;
plot( f_CPC, P_CPC ); hold on;
plot( f_IPC, P_IPC ); grid on; hold off;
xlabel( 'Frequency, $\left( Hz \right)$' )
ylabel( 'Amplitude, $\left( N.m \right)$' )

xlim( [0.00 0.25] );
ylim( [0.00 10e+14] );

title( 'PSD of BRBM' );

% Add legend
legend( 'CPC', 'CPC+IPC', ...
        'NumColumns'    , 1             , ...
        'Location'      , 'northeast' );
make_nice_plot( PRNT, figs_dir, 'PSD' );

% Get percent difference
max_psd_CPC = max( P_CPC(2:end) );
max_psd_IPC = max( P_IPC(2:end) );
prct_diff_psd = abs(max_psd_CPC - max_psd_IPC)/max_psd_CPC


%% Plot

labels      = [ {'Rotor speed'      , '$\left( rad/s \right)$'}; ...
                {'CP Angle'         , '$\left( deg   \right)$'}; ...
                {'BRBM'             , '$\left( N.m   \right)$'}; ...
                {'IP Angle'         , '$\left( deg   \right)$'}; ...
                {'Generator Torque' , '$\left( N.m   \right)$'}; ...
                {'Turbine Power'    , '$\left( W     \right)$'} ];
ROWS    = 3; COLS = 2;               	% Subplot number of desired rows and columns
IDX     = 1:ROWS*COLS;               	% Index for subplots

figure( CNTR ); CNTR = CNTR + 1;
for i=1:6
    % Create subplot space
    subplot( ROWS, COLS, IDX(i) )

    % Plot entries
    plot( CPC_Data{1}, CPC_Data{i+1} ); grid on; hold on;
    if( i==2 || i==4 )
        yyaxis right;
        plot( IPC_Data{1}, IPC_Data{i+1} ); hold off;
        % Add labels
        yyaxis left;
        xlabel( 'Time $\left( s     \right)$' );
        ylabel( [labels(i,1), labels{i,2}] );
    else
        plot( IPC_Data{1}, IPC_Data{i+1} ); hold off;
        % Add labels
        xlabel( 'Time $\left( s     \right)$' );
        ylabel( [labels(i,1), labels{i,2}] );
    end

    % Correct the x-lim
    xlim( [0 750] );
    
    % Add legend on only the first subplot
    if( i == 1 )
        legend( 'CPC', 'CPC+IPC', ...
                'NumColumns'    , 1         , ...
                'Location'      , 'southeast' );
    end

    % Conditional formatting
    if( i == 1 )
        ylim( [0.77 0.80] );
        yticks( [0.77 0.78 0.79 0.80] )
        yticklabels( [ "0.77" "0.78" "0.79" "0.80" ] );

    elseif( i == 2 )
        yyaxis left;
        ylim( [0 6] );
        yticks( [0 2 4 6] )
        yticklabels( [ "0" "2" "4" "6" ] );

        yyaxis right;
        ylim( [0.2 0.8] );
        yticks( [0.2 0.4 0.6 0.8] )
        yticklabels( [ "0.20" "0.40" "0.60" "0.80" ] );

    elseif( i == 3 )
        ylim( [4 7].*1e7 );
        yticks( [4 5 6 7].*1e7 )
        yticklabels( [ "4.00" "5.00" "6.00" "7.00" ] );
        text(0, 1.00    ,	'$\times 10^{7}$'  	, ...
                    'horiz'     ,	'left'              , ...
                    'vert'      ,	'bottom'            , ...
                    'Units'     ,	'normalized'        , ...
                    'fontsize'  ,	fontsize-1 );

    elseif( i == 4 )
        yyaxis right;
        ylim( [1.90 1.96] );
        yticks( [ 1.90 1.92 1.94 1.96] )
        yticklabels( [ "1.90" "1.92" "1.94" "1.96" ] );

    elseif( i == 5 )
        ylim( [1.86 1.92].*1e7 );
        yticks( [1.86 1.88 1.90 1.92].*1e7 )
        yticklabels( [ "1.86" "1.88" "1.90" "1.92" ] );
        text(0, 1.00    ,	'$\times 10^{7}$'  	, ...
                    'horiz'     ,	'left'              , ...
                    'vert'      ,	'bottom'            , ...
                    'Units'     ,	'normalized'        , ...
                    'fontsize'  ,	fontsize-1 );

    elseif( i == 6 )
        ylim( [12 18].*1e6 );
        yticks( [12 14 16 18].*1e6 )
        yticklabels( [ "12" "14" "16" "18" ] );
        text(0, 1.00    ,	'$\times 10^{6}$'  	, ...
                    'horiz'     ,	'left'              , ...
                    'vert'      ,	'bottom'            , ...
                    'Units'     ,	'normalized'        , ...
                    'fontsize'  ,	fontsize-1 );
    end

end
sgtitle( 'CPC vs CPC+IPC' );
make_nice_plot( PRNT, figs_dir, 'comparison' );


%% Store updated mat file
% 
% storedStructure = load( fullfile( data_dir, 'QFT_AzimuthVar_4x4IndepSISO_Simulink_HI_GAIN.mat') );
% Aclass = storedStructure.Aclass;
% names = storedStructure.names;
% data = storedStructure.data;
% 
% A = max( abs(data(:,3)) );
% data(:,3) = data(:,3)*1/2 + A*0.470;
% 
% mat_name = fullfile( data_dir, 'QFT_AzimuthVar_4x4IndepSISO_Simulink_LOWER_GAIN.mat');
% save( mat_name, 'Aclass', 'names', 'data');

%% Load data

