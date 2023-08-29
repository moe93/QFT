% Try out different responses for specification 7
%   Reference tracking specification
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

%% Specification 7 (Type 6 in Dr. Garcia's book)
% --- Type 6: Reference tracking specification
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
%   ***NOTE ON (eps_U) AND (eps_L):
%           eps_U and eps_L control the steady-state error. For instance, a
%       value of (eps_U = 0.01) will result in a steady-state response of
%       (y_U = 1.01) when subjected to a step input.
%           Similarily, a value of (eps_L = 0.15) will result in a
%       steady-state response of (y_L = 0.85) when subjected to a step
%       input.
%
%   ***NOTE ON (a_U) AND (a_L):
%           The value of the pole (a_*) influences the rise time and
%       settling time of the system response. The higher the (a_*) value,
%       the lower the the rise time and settling time of the system
%       response.
%           In other words, if you want the system to have a fast reponse,
%       choose a high value for the pole (a_*) value. Similarly, if you
%       want to have a slower system response, choose a low value for the
%       pole (a_*).
%

% Frequencies of interest
% omega_6 = [ 1e-2 1e-1 1e0 1e1 ];
omega_6 = [ 1e-1 1e0 ];

% Restriction
% Upper bound
% a_U = 0.15; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.01;
a_U = 0.10; zeta = 0.8; wn = 1.25*a_U/zeta; eps_U = 0.01;
num = [ conv([1/a_U 1], [0 1+eps_U]) ];
den = [ (1/wn)^2 (2*zeta/wn) 1 ];
del_6_hi = tf( num, den );

% Lower bound
% a_L = 0.25; eps_L = 0.00;
a_L = 0.35; eps_L = 0.01;
num = 1-eps_L;
den = [ conv([1/a_L 1], [1/a_L 1]) ];
del_6_lo = tf( num, den );

% Tracking weight
del_6 = [ del_6_hi  ;
          del_6_lo ];

% --- Plot bounds
figure( CNTR ); CNTR = CNTR + 1;
step( del_6(1) );   hold on ;  grid on;
step( del_6(2) );   hold off;
legend( "Upper Bound", "Lower Bound" )
make_nice_plot();
