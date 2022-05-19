function strak = strakImport( filename )
% fng_var14_33stn.hsf
% strak = strakImport( 'fng_var14_33stn.hsf' )

str = fileread( filename );

idx_line_end = strfind( str, newline );
num_lines = length( idx_line_end );

strak.eta = [];
strak.x = {};
strak.z = {};

idx_current_eta = 0;
current_start_line_idx = 0;
idx_of_current_eta = 0;

for current_line_idx = 1:num_lines
    
    if current_line_idx == 1
        current_idx = 1;
    else
        current_idx = idx_line_end(current_line_idx-1)+1;
    end
    current_line_content = str(current_idx:idx_line_end(current_line_idx));
    
    if current_line_content(1) == 'S'
        idx_current_eta = idx_current_eta + 1;
        idx_eta_in_line = strfind( current_line_content, 'eta = ' );
        idx_semicolon_in_line = strfind( current_line_content, ';' );
        strak.eta(idx_current_eta) = str2double( current_line_content( (idx_eta_in_line+6):(idx_semicolon_in_line-1) ) );
        current_start_line_idx = current_line_idx + 3;
        idx_of_current_eta = 1;
    end
    
    if current_line_idx == current_start_line_idx - 2
        if length(strak.eta) == 1
            span_in = str2double( current_line_content( 2:10 ) );
        else
            span_tmp = str2double( current_line_content( 2:10 ) );
        end
    end
    
    if current_line_idx >= current_start_line_idx && length(current_line_content) >= 20
        if length(strak.x) < idx_current_eta
            strak.x{ idx_current_eta } = str2double( current_line_content(2:10) );
            strak.z{ idx_current_eta } = str2double( current_line_content(13:20) );
        else
            strak.x{ idx_current_eta }(end+1) = str2double( current_line_content(2:10) );
            strak.z{ idx_current_eta }(end+1) = str2double( current_line_content(13:20) );
        end
        idx_of_current_eta = idx_of_current_eta + 1;
    end
    
end

strak.span = ( span_tmp - span_in ) / ( strak.eta(end) - strak.eta(1) );

end