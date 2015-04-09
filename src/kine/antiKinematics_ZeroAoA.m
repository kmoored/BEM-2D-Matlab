function [ xp_new,zp_new ] = antiKinematics_ZeroAoA( xp,zp,h_c,c,f,t,Qinf,mesh_type,ramped,i_t )
    Npanels = length(xp)-1;
    Nx = Npanels/2 + 1;

    if strcmp(mesh_type,'fluid')
        xpt = xp(Nx:end);
    elseif strcmp(mesh_type,'solid')
        xpt = xp(1:end);
    else
        fprintf('ERROR! Incorrect mesh type ''%s.'' Valid mesh types are ''fluid'' or ''solid.''',mesh_type)
        return
    end

    % Rotate back to zero degree attack
    theta_rot = atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf);
    if strcmp(mesh_type,'fluid')
        xrot(1:Npanels+1,1) = xp(0.5*Npanels+1,1);
        zrot(1:Npanels+1,1) = zp(0.5*Npanels+1,1);
    elseif strcmp(mesh_type,'solid')
        xrot(1:Npanels+1,1) = xp(1,1);
        zrot(1:Npanels+1,1) = zp(1,1);
    else
        fprintf('ERROR! Incorrect mesh type ''%s.'' Valid mesh types are ''fluid'' or ''solid.''',mesh_type)
        return
    end
    xp_new = ((xp-xrot).*cos(-1*theta_rot) - (zp-zrot).*sin(-1*theta_rot)) + xrot;
    zp_new = ((zp-zrot).*cos(-1*theta_rot) + (xp-xrot).*sin(-1*theta_rot)) + zrot;

    % Translate the heaved LE position back to the origin
    if strcmp(mesh_type,'fluid')
        xp_new = xp_new - xp_new(0.5*Npanels+1,1);
        zp_new = zp_new - ramped(i_t)*(h_c*c*sin(2*pi*f*t));
    elseif strcmp(mesh_type,'solid')
        xp_new = xp_new - xp_new(1,1);
        zp_new = zp_new - ramped(i_t)*(h_c*c*sin(2*pi*f*t));
    else
        fprintf('ERROR! Incorrect mesh type ''%s.'' Valid mesh types are ''fluid'' or ''solid.''',mesh_type)
        return
    end


end

