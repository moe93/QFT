function print_f( f_ii )
%PRINT_G Summary of this function goes here
%   Detailed explanation goes here

f = f_ii;

f_num = cell2mat( tf(f).num );
fprintf( "NUMERATOR\n" )
for i=1:length( f_num )
    fn = f_num(i);
    if( i == 1 )
        % --- Conditional formatting
        if( (abs(fn) >= 1e4 || abs(fn) <= 1e-4) )
            fprintf( "{ %3.3e, ", fn )
        else
            fprintf( "{ %3.4f, ", fn )
        end

    elseif( i > 1 && i ~= length(f_num) )
        % --- Conditional formatting
        if( abs(fn) >= 1e4 || abs(fn) <= 1e-4 )
            fprintf( "%3.3e, ", fn )
        else
            fprintf( "%3.4f, ", fn )
        end

    else
        % --- Conditional formatting
        if( abs(fn) >= 1e4 || abs(fn) <= 1e-4 )
            fprintf( "%3.3e }\n\n", fn )
        else
            fprintf( "%3.4f }\n\n", fn )
        end
    end
end

f_den = cell2mat( tf(f).den );
fprintf( "DENOMERATOR\n" )
for i=1:length( f_den )
    fd = f_den(i);
    if( i == 1 )
        % --- Conditional formatting
        if( (abs(fd) >= 1e4 || abs(fd) <= 1e-4) )
            fprintf( "{ %3.3e, ", fd )
        else
            fprintf( "{ %3.4f, ", fd )
        end

    elseif( i > 1 && i ~= length(f_den) )
        % --- Conditional formatting
        if( abs(fd) >= 1e4 || abs(fd) <= 1e-4 )
            fprintf( "%3.3e, ", fd )
        else
            fprintf( "%3.4f, ", fd )
        end

    else
        % --- Conditional formatting
        if( abs(fd) >= 1e4 || abs(fd) <= 1e-4 )
            fprintf( "%3.3e }\n\n", fd )
        else
            fprintf( "%3.4f }\n\n", fd )
        end
    end
end