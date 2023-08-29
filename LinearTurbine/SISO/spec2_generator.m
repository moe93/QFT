% Try out different responses for specification 2
%   Sensitivity or ouptut disturbance rejection specification
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Jul. 28th, 2023
%
% CHANGELOG :
%   Jul. 28th, 2023
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

%% Specification 2 (Type 3 in Dr. Garcia's book)
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
% omega_3 = [ 1e-2 1e-1 1e0 1e1 ];
omega_3 = [ 1e-1 1e0 ];

% Restriction
a_d     = 0.01;
num     = [ 1/a_d   , 0 ];
den     = [ 1/a_d   , 1 ];
% num     = [ 0.025   , 0.2   , 0.018 ];
% den     = [ 0.025   , 10    , 1     ];
del_3   = tf( num, den );

% --- Plot bounds
figure( CNTR ); CNTR = CNTR + 1;            % Instantiate handle
w_3 = logspace( log10(min(omega_3)), ...    % Generate freq. vector
                log10(max(omega_3)), ...    % ...
                1024 );                     % ...
bode( del_3, w_3 ); grid on;                % Plot Bode diagram
make_nice_plot();
