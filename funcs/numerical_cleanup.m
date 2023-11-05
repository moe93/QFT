function P = numerical_cleanup( P, threshold, TOL )
% Remove value below a certain threshold to reduce numerical instability
%   For example, values below 1e-12 can be safely replaced with 0.

arguments
    P           (:,:,:) {mustBeA(P,'tf')} = false
    threshold   (1,1)   {mustBeNumeric}   = 1e-08
    TOL         (1,1)   {mustBeNumeric}   = 0.01
end

% if( nargin >= 1 && nargin < 3 )
%     error( "If PRNT == true, a directory AND name must be provided" );
% end

% --- Get total plants size
% [x-dim, y-dim, z-dim] = [nrowsP, ncolsP, nvarsP]
[~, ~, nvarsP] = size( P );

% --- Cleanup plants transfer function by removing values below 1e-08
for ii = 1:nvarsP
    [n, d] = tfdata( minreal(P( :, :, ii, 1 ), TOL) );
    n = cellfun(@(x) {x.*(abs(x) > threshold)}, n);
    d = cellfun(@(x) {x.*(abs(x) > threshold)}, d);
    P( :, :, ii, 1 ) = tf(n, d);
end

end