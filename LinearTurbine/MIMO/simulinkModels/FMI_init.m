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

t_linearization = [ 163.1450  163.8150  164.4850  165.1550  165.8250  ...
                    166.4950  167.1650  167.8350  168.5050  169.1750  ...
                    169.8440  170.5140  171.1830 ];

% open_system( '.\Simulink_Linearization_Model' );
open_system( '.\Simulink_Linearization_Model_FMIKit' );

