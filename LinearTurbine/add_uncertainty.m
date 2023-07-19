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
inputNames  = [ "u_{pitch}" ];
outputNames = [ "omega" ];
% State-space model
sys         = ss( A, B, C, D                , ...
                  'StateName' , stateNames  , ...
                  'InputName' , inputNames  , ...
                  'OutputName', outputNames );

% --- Generate TF from SS model
TF = tf( sys );

%% Manaully construct

A2 = [ 0 1; -3.123e-8 -0.076995 ];
B2 = [ 0 -1.535 ].';
C2 = [ 0 1 ];
D2 = 0;

sys2         = ss( A2, B2, C2, D2           , ...
                  'StateName' , stateNames  , ...
                  'InputName' , inputNames  , ...
                  'OutputName', outputNames );
% --- Generate TF from SS model
TF2 = tf( sys2 );


%% Manaully construct SISO

A3 = [ -0.076995 ];
B3 = [ -1.535 ].';
C3 = [ 1 ];
D3 = 0;

sys3         = ss( A3, B3, C3, D3 );
% --- Generate TF from SS model
TF3 = tf( sys3 );

%% Check plant against Nyquist stability guidelines

P0 = TF3;
output = nyquistStability( P0 );

if( PLOT )
    figure();    bode( P0 ); grid on;
    figure();  rlocus( P0 ); grid on;
    figure(); nichols( P0 ); grid on;
    figure(); nyquist( P0 );
end

%% Check influence of uncertainty on response
% Generate SS model using analytical approach

% M_0 = 0.5   ;                   % Mass of cart                  [  kg  ]
% m_0 = 0.2   ;                   % Mass of rod                   [  kg  ]
% b_0 = 0.1   ;                   % Co-efficient of friction
% I_0 = 0.006 ;                   % 2nd mass moment of inertia    [kg.m^2]
% g_0 = 9.80665 ;                 % Gravitational acceleration    [m.s^-2]
% l_0 = 0.6/2 ;                   % Half the rod length           [  m   ]
% 
% p = I_0*(M_0+m_0)+M_0*m_0*l_0^2;    % Denominator for the A and B matrices
% 
% A_theory = [ 0              1                       0               0   ;
%              0  -(I_0+m_0*l_0^2)*b_0/p  (m_0^2*g_0*l_0^2)/p         0   ;
%              0              0                       0               1   ;
%              0  -(m_0*l_0*b_0)/p        m_0*g_0*l_0*(M_0+m_0)/p     0 ] ;
% B_theory = [        0            ;
%              (I_0+m_0*l_0^2)/p   ;
%                     0            ;
%                 m_0*l_0/p       ];
% C_theory = [ 1  0   0   0  ;
%              0  0  -1   0 ];
% D_theory = [ 0  ;
%              0 ];
% 
% % States and inputs names
% stateNames  = [ "phi" "omega" ];
% inputNames  = [ "u_{pitch}" ];
% outputNames = [ "omega" ];
% % State-space model
% sys_theory = ss( A_theory, B_theory, C_theory, D_theory , ...
%                  'StateName'    ,   stateNames          , ...
%                  'InputName'    ,   inputNames          , ...
%                  'OutputName'   ,   outputNames );
% 
% % --- Generate TF from SS model
% TF_theory = tf( sys_theory );
