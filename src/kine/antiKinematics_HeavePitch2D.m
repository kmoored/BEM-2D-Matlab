function [xp_new,zp_new] = antiKinematics_HeavePitch2D(xp,zp,alpha_max,h_c,f,t,phi,mesh_type)
    Npanels = length(xp)-1;
    Nx = Npanels/2 + 1;

    if strcmp(mesh_type,'fluid')
        xpt = xp(Nx:end);
        %c = xpt(end) - xpt(1);
        c = 1;
    elseif strcmp(mesh_type,'solid')
        xpt = xp(1:end);
        %c = xpt(end) - xpt(1);
        c = 1;
    else
        fprintf('ERROR! Incorrect mesh type ''%s.'' Valid mesh types are ''fluid'' or ''solid.''',mesh_type)
        return
    end

    % Rotate back to zero degree attack
    theta_rot = alpha_max*sin(2*pi*f*t+phi);
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
        zp_new = zp_new - h_c*c*sin(2*pi*f*t);
    elseif strcmp(mesh_type,'solid')
        xp_new = xp_new - xp_new(1,1);
        zp_new = zp_new - h_c*c*sin(2*pi*f*t);
    else
        fprintf('ERROR! Incorrect mesh type ''%s.'' Valid mesh types are ''fluid'' or ''solid.''',mesh_type)
        return
    end
end

