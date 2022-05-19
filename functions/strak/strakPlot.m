function strakPlot( strak, varargin )

% strakPlot(strak,20000)

if isempty(varargin)
    is_plot_chord = false;
else
    is_plot_chord = varargin{1};
end

num_eta = length(strak.eta);

for i = 1:num_eta
    
    len_section = length( strak.x{i} );
    
    plot3( strak.x{i}, repmat(strak.span*strak.eta(i),1,len_section), strak.z{i} )
    hold on
    
    if isfield(strak,'xz_lead') && is_plot_chord
%         xyz_trail = [ strak.xz_trail{i}(1), strak.span*strak.eta(i), strak.xz_trail{i}(2) ];
%         xyz_lead = [ strak.xz_lead{i}(1), strak.span*strak.eta(i), strak.xz_lead{i}(2) ];
%         plot3( [xyz_lead(1),xyz_trail(1)], [xyz_lead(2),xyz_trail(2)], [xyz_lead(3),xyz_trail(3)], 'k-' )
        len_c = length( strak.camber_line_xz{i}(1,:) );
        plot3( strak.camber_line_xz{i}(1,:), repmat( strak.span*strak.eta(i), 1, len_c ), strak.camber_line_xz{i}(2,:), 'k--' )
    end
    
end

axis equal

end