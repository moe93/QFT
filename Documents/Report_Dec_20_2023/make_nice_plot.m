function make_nice_plot( PRNT, dir, name )
% Make the plot look cleaner and nicer
% in addition, prepare to print to PDF

arguments
    PRNT    (1,1) {mustBeNumericOrLogical}  = false
    dir     (1,:) string  {mustBeText}      = 'foo'
    name    (1,1) string  {mustBeText}      = 'bar'
end

if( nargin >= 1 && nargin < 3 )
    error( "If PRNT == true, a directory AND name must be provided" );
end

% Reformat to look nice and stuff
set( gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4' )
screen   = get( 0   , 'ScreenSize' );
A4       = get( gcf , 'PaperSize'  );
ratio    = A4(2)/A4(1);
position = [ 60             (screen(4)-nearest((screen(3)-600)/ratio)-120) ...
             screen(3)-600  nearest((screen(3)-600)/ratio) ];
set( gcf, 'Units'   , 'pixels' );
set( gcf, 'Position', position );

paper_offset = 0.5; % centimeters

% Prepare for printing to a PDF
set( gcf, 'PaperType', 'A4', 'PaperOrientation', 'Landscape' )
paper       = get( gcf, 'PaperSize' );
paper_pos   = [ paper_offset                paper_offset            ...
                paper(1)-2*paper_offset     paper(2)-2*paper_offset ];
set(gcf,'PaperPosition',paper_pos)

if( PRNT )                                      % If flag is set
    if( ~isfolder( dir ) ); mkdir( dir ); end   % Check if directory exists
    
    name = fullfile( dir, name );
    % print( gcf, '-dsvg', name, '-r600' )        %   Print .svg file
    print( gcf, '-dpng', name, '-r600' )	    %   Print .png  file
    close;                                      %   Close windoq
else                                            % Else
    0;                                          %   Do nothing
end

end

