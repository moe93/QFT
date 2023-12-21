function [f, P1] = do_fft( data, time, sampling_time, x_lim, y_lim )
%do_fft( data, time, sampling_time )
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
Y = fft( x );
P2 = abs(Y/L);

if( mod(L, 2) ~= 0 )                            % If NOT even
    P1 = P2(1:floor(L/2+1));                    %   Floor the index
else                                            % Else
    P1 = P2(1:L/2+1);                           %   Do not manipulate
end
P1(2:end-1) = 2*P1(2:end-1);                    % Single-Sided Amplitude Spectrum
f = Fs*(0:(L/2))/L;

if( ~nargout )
    % Plot
    subplot( 2, 1, 1 );
    plot( t, x, 'color', [0.0000 0.4470 0.7410] ); hold on; grid on;
    
    subplot( 2, 1, 2 );
    plot( f, P1, 'color', [0.8500 0.3250 0.0980] ); hold on; grid on;
    xlim( x_lim );
    
    xlabel( 'Frequency $\left( Hz \right)$', ...  	% x-axis label
            'interpreter', 'latex' );               % ...
    ylabel( 'Amplitude', ...                        % y-axis label
            'interpreter', 'latex' );               % ...
    txt = ['Single-Sided Amplitude Spectrum'];   	% Title text
    title( txt, 'interpreter', 'latex' );       	% Add title
    % hold off;
end

end

