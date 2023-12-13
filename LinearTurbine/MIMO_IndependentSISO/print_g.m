function print_g( g_ii )
%PRINT_G Summary of this function goes here
%   Detailed explanation goes here

g = g_ii;

g_num = cell2mat( tf(g).num );
fprintf( "NUMERATOR\n" )
for i=1:length( g_num )
    gn = g_num(i);
    if( i == 1 )
        % --- Conditional formatting
        if( (abs(gn) >= 1e4 || abs(gn) <= 1e-4) )
            fprintf( "{ %3.3e, ", gn )
        else
            fprintf( "{ %3.4f, ", gn )
        end

    elseif( i > 1 && i ~= length(g_num) )
        % --- Conditional formatting
        if( abs(gn) >= 1e4 || abs(gn) <= 1e-4 )
            fprintf( "%3.3e, ", gn )
        else
            fprintf( "%3.4f, ", gn )
        end

    else
        % --- Conditional formatting
        if( abs(gn) >= 1e4 || abs(gn) <= 1e-4 )
            fprintf( "%3.3e }\n\n", gn )
        else
            fprintf( "%3.4f }\n\n", gn )
        end
    end
end

g_den = cell2mat( tf(g).den );
fprintf( "DENOMERATOR\n" )
for i=1:length( g_den )
    gd = g_den(i);
    if( i == 1 )
        % --- Conditional formatting
        if( (abs(gd) >= 1e4 || abs(gd) <= 1e-4) )
            fprintf( "{ %3.3e, ", gd )
        else
            fprintf( "{ %3.4f, ", gd )
        end

    elseif( i > 1 && i ~= length(g_den) )
        % --- Conditional formatting
        if( abs(gd) >= 1e4 || abs(gd) <= 1e-4 )
            fprintf( "%3.3e, ", gd )
        else
            fprintf( "%3.4f, ", gd )
        end

    else
        % --- Conditional formatting
        if( abs(gd) >= 1e4 || abs(gd) <= 1e-4 )
            fprintf( "%3.3e }\n\n", gd )
        else
            fprintf( "%3.4f }\n\n", gd )
        end
    end
end