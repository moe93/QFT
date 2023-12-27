function [f, psdx] = do_psd( data, time, sampling_time, x_lim, y_lim )
%do_psdx( data, time, sampling_time )
%   Conactenate the FFT generation process into a custom function
%   INPUTS:-
%       - data          : Data vector (row OR column) we wish to perform FFT on
%       - time          : Data corresponding time vector (row OR column)
%       - sampling_time : Scalar corresponding to the sampling time of the
%                           data vector
%
%   RETURN:-
%       - f     : Freuquency content found in data
%       - P1    : Single-Sided Amplitude Spectrum
%
arguments
    data            ( :, 1 )    double
    time            ( :, 1 )    double
    sampling_time               double
    x_lim                       double = [0, 1]
    y_lim                       double = [0, 1]
end

% Unpack data
x = data;                                       % Data to run through FFT
t = time;                                       % Data corresponding time vector
Ts = sampling_time;                             % Data sampling time

% Perform FFT 
L = length( data );                          	% Length of data set
Fs = 1/Ts;                                      % Sampling frequency
xdft = fft( x );

if( mod(L, 2) ~= 0 )                            % If NOT even
    xdft = xdft(1:floor(L/2+1));                %   Floor the index
else                                            % Else
    xdft = xdft(1:L/2+1);                       %   Do not manipulate
end
psdx = (1/(Fs*L)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
f = 0:Fs/length(x):Fs/2;

if( ~nargout )
    % Plot
    subplot( 2, 1, 1 );
    plot( t, x, 'color', [0.0000 0.4470 0.7410] ); hold on; grid on;
    
    subplot( 2, 1, 2 );
    plot( f, psdx, 'color', [0.8500 0.3250 0.0980] ); hold on; grid on;
    xlim( x_lim ); ylim( y_lim )
    
    xlabel( 'Frequency $\left( Hz \right)$', ...  	% x-axis label
            'interpreter', 'latex' );               % ...
    ylabel( 'Amplitude', ...                        % y-axis label
            'interpreter', 'latex' );               % ...
    txt = ['Power Spectral Density (PSD)'];   	    % Title text
    title( txt, 'interpreter', 'latex' );       	% Add title
    % hold off;
end

end

