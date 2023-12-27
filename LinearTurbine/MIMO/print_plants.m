function print_plants( P )
%PRINT_PLANTS Format print plants for easier perusing
%   This function's main purpose is NOT to give accurate representation of
%   the plants' transfer functions (TFs), but rather give a simple idea of
%   the plants structure.

[xP, yP, zP] = size( P );

syms s;
for z = 1:zP
    fprintf( "Variation (%i)\n", z );
    fprintf( "==============\n\n" );
    for x = 1:xP
        fprintf( '\tP( %i, %i )\n', x, x );
        P_ii = P(x, x, z); 
        num_txt = string( vpa(poly2sym(cell2mat(P_ii.num), s), 3) );
        den_txt = string( vpa(poly2sym(cell2mat(P_ii.den), s), 3) );

        div_len = max( strlength(num_txt), strlength(den_txt));
        div_txt = repmat( '-', 1, div_len );

        fprintf( "\t\t%s\n"  , num_txt );
        fprintf( "\t\t%s\n"  , div_txt );
        fprintf( "\t\t%s\n\n", den_txt );
    end
    fprintf( '\n' );
end

end


% for j=1:length(a)
%     fprintf( "\t\tNat'l freq.  w_%i = %0.3f (rad/s) OR %0.3f (Hz).", ...
%                 j, omgr(j), omgr(j)/(2*pi) );
%     text =   "\t\tMode shape {X}_%i = {%s}.\n";
%     text = sprintf( text, j, repmat('%0.3f ', 1, length(a)) );
%     text = sprintf( text, a(:,j) ); disp( text );
% end
