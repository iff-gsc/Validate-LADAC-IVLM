function strak = strakSetCamber( strak )

num_section = length(strak.eta);

num_camber_points = 100;

for i = 1:num_section
    p = polyfit(strak.x{i},strak.z{i},1);
    alpha_regression = atan( -p(1) );
    xz_rot = rotationMatrix( alpha_regression ) * [ strak.x{i}; strak.z{i} ];
    [~,idx_max] = max( xz_rot(1,:) );
    [~,idx_min] = min( xz_rot(1,:) );
    strak.xz_lead{i} = [strak.x{i}(idx_min);strak.z{i}(idx_min)];
    strak.xz_trail{i} = [strak.x{i}(idx_max);strak.z{i}(idx_max)];
    strak.alpha(i) = atan( (strak.xz_lead{i}(2)-strak.xz_trail{i}(2))/-(strak.xz_lead{i}(1)-strak.xz_trail{i}(1)) );
    xz_rot = rotationMatrix( strak.alpha(i) ) * [ strak.x{i}; strak.z{i} ];
    chord_length = xz_rot(1,idx_max) - xz_rot(1,idx_min);
    if idx_min < idx_max
        idx_1 = idx_min:idx_max;
        idx_2 = [idx_min:-1:1,length(xz_rot(1,:)):-1:idx_max];
    else
        idx_1 = [idx_min:length(xz_rot(1,:)),1:idx_max];
        idx_2 = idx_min:-1:idx_max;
    end
    strak.camber_line_xz{i} = zeros(2,num_camber_points);
    strak.camber_line_xz{i}(1,:) = xz_rot(1,idx_1(1)):chord_length/(num_camber_points-1):xz_rot(1,idx_1(end));
    z_1 = interp1( xz_rot(1,idx_1), xz_rot(2,idx_1), strak.camber_line_xz{i}(1,:) );
    z_2 = interp1( xz_rot(1,idx_2), xz_rot(2,idx_2), strak.camber_line_xz{i}(1,:) );
    strak.camber_line_xz{i}(2,:) = (z_1+z_2)/2;
    z_mean = strak.camber_line_xz{i}(2,:) - xz_rot(2,1);
    [strak.camber(i),idx_camber_max] = max( abs( z_mean ) );
    strak.camber(i) = strak.camber(i) * sign( z_mean(idx_camber_max) ) / abs( chord_length );
    strak.Xf(i) = ( strak.camber_line_xz{i}(1,idx_camber_max) - strak.camber_line_xz{i}(1,1) ) / chord_length;
    if strak.Xf(i) < 0
        strak.Xf(i) = 1 + strak.Xf(i);
    end
    if abs( xz_rot(2,idx_min) - xz_rot(2,idx_max) ) > 1e-8
        error('transformation did not work correctly')
    end
    strak.camber_line_xz{i} = rotationMatrix( -strak.alpha(i) ) * strak.camber_line_xz{i};
end

end

function rot_mat = rotationMatrix( alpha )
    rot_mat = [ cos(alpha), -sin(alpha); sin(alpha), cos(alpha) ];
end