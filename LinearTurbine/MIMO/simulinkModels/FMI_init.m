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

% t_linearization = [ 163.1450  163.8150  164.4850  165.1550  165.8250  ...
%                     166.4950  167.1650  167.8350  168.5050  169.1750  ...
%                     169.8440  170.5140  171.1830 ];
t_linearization = [ 188.3960  188.9640  189.5320  190.0980  190.6660  ...
                    191.2340  191.8020  192.3700  192.9360  193.5040  ...
                    194.0720  194.6400  195.2060 ];

%% Simple results analysis

idx = zeros( length(t_linearization), 1 );
for i = 1:length( t_linearization )
    t = t_linearization( i );
    idx( i ) = find( abs(out.Azimuth_angle.Time - t) <= 1e-3, 1, 'first' );
end

out.Azimuth_angle.Data( idx )

%% Save results to .mat file
%

save Simulink_Linearized_Model_REV002.mat out Simulink_Linearization_Model_FMIKit_Timed_Based_Linearization;
