% FMIKit Loader
% 
% Author: Mohammad Odeh
% Date  : Dec.  1st, 2023

%% Setup environment

clc;    % Clear screen

% --- Controllers and pre-filter directories
scriptDir       = mfilename( 'fullpath' );
scriptPathParts = strsplit( scriptDir, filesep );

% Name of directory where FMIKit is located (default folder name)
FMI_ver     = '3.1';
FMI_name    = ['FMIKit-Simulink-' FMI_ver];
% Go up one level and generate new path
FMI_dir     = fullfile( scriptPathParts{1:end-4}, FMI_name );

% Check if directory exists, if not, download and unzip from GitHub
if ~isfolder( FMI_dir )
    fprintf( '[INFO] FMIKit NOT found. Downloading %s from GitHub...', FMI_name );
    FMI_url = ['https://github.com/CATIA-Systems/FMIKit-Simulink/' ...
               'releases/download/v' FMI_ver '/' FMI_name '.zip'];
    unzip( FMI_url, FMI_name );
    fprintf( 'DONE!\n' );
else
    fprintf( '[INFO] %s found.\n', FMI_name );
end

% Add FMIKit to path
addpath( FMI_dir );

% Initialize
FMIKit.initialize()

% Clean-up workspace
clear FMI* script*

%% Load Model
%

% open_system( '.\Simulink_Linearization_Model' );
open_system( '.\Simulink_Linearization_Model_FMIKit' );

%% Required model parameters
%

% % [0:30:360]
% t_linearization = [ 188.3960  188.9640  189.5320  190.0980  190.6660  ...
%                     191.2340  191.8020  192.3700  192.9360  193.5040  ...
%                     194.0720  194.6400  195.2060 ];

% % [0:15:360]
% t_linearization = [ 188.3960  188.6800  188.9640  189.2470  189.5320  ...
%                     189.8140  190.0980  190.3800  190.6660  190.9500  ...
%                     191.2340  191.5180  191.8020  192.0860  192.3700  ...
%                     192.6520  192.9360  193.2200  193.5040  193.7880  ...
%                     194.0720  194.3560  194.6400  194.9240  195.2060 ].';

% [0:15:360)
t_linearization = [ 188.3960  188.6800  188.9640  189.2470  189.5320  ...
                    189.8140  190.0980  190.3800  190.6660  190.9500  ...
                    191.2340  191.5180  191.8020  192.0860  192.3700  ...
                    192.6520  192.9360  193.2200  193.5040  193.7880  ...
                    194.0720  194.3560  194.6400  194.9240 ].';

%% Simple results analysis

idx = zeros( length(t_linearization), 1 );
for i = 1:length( t_linearization )
    t = t_linearization( i );
    idx( i ) = find( abs(out.Azimuth_angle.Time - t) <= 1e-3, 1, 'first' );
end

angles_linearization = out.Azimuth_angle.Data( idx );
angles = table( t_linearization, round(angles_linearization) )

%% Save results to .mat file
%

save Simulink_Linearized_Model_REV004.mat out Simulink_Linearization_Model_FMIKit_Timed_Based_Linearization;
