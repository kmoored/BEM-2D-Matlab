

function [xstar_n,zstar_n] = fencing(xpts,zpts,up,wp,xp,zp,vn,x_b,z_b,delT,epBod)

    N = 10;
%     Nbod = 1000;
    
    xpts = xpts - x_b;
    zpts = zpts - z_b;
    
    xp = xp - x_b;
    zp = zp - z_b;
    
    L = sqrt((xp(2:end) - xp(1:end-1)).^2 + (zp(2:end) - zp(1:end-1)).^2);
    Lmax = max(L);
    
    xstar_n = xpts;
    zstar_n = zpts;

  
    % Defining the top and bottom surfaces
    Npan = length(xp) - 1;
%     xb = xp(1:Npan/2 + 1);
%     xt = xp(Npan/2 + 1:end);
%     
%     zb = zp(1:Npan/2 + 1);
%     zt = zp(Npan/2 + 1:end);
    
    % Defining the outward normal vectors at each panel endpoint.
    vnp = zeros(length(xp),2);
    vnp(1,:) = vn(1,:)/2 + vn(end,:)/2;
    vnp(end,:) = vn(1,:)/2 + vn(end,:)/2;
    vnp(2:end-1,:) = vn(1:end-1,:)/2 + vn(2:end,:)/2;
    


    for i = 1:length(xpts)
        
        % Extracting panel starting coordinates
        x0 = xpts(i);
        z0 = zpts(i);
        
        % Extracting velocity field at starting coordinates
        u = up(i);
        w = wp(i);
        
        % Calculating the advection of the panel points (x0,z0)-->(x',z')
        xpri = x0 + u*delT;
        zpri = z0 + w*delT;
        
        % Calculating the tangential and normal vectors to the direction
        % [u,w]'.
        tvec = [xpri - x0 zpri - z0];
        t = tvec/norm(tvec);
        n = [-t(2) t(1)];
        R0 = [t; n];
        
        % Calculating the vector from the point (x',z') to all of the
        % body panel points. Calculating the distance and the unit vector
        % r_hat.
        r = [xp - xpri zp - zpri];
        D = sqrt(r(:,1).^2 + r(:,2).^2);
        r_hat = diag(1./D)*r;
        
        
        % Finding minimum distance point (xmin,zmin)
        [~,indmin] = min(D);
        
        % Calculating (n_hat dot r_hat)
        Val = dot(vnp(indmin,:),r_hat(indmin,:));
        
        % Determining whether the point (x',z') is inside or outside of the
        % body.
        if Val > 0
            fence = 1;
        else
            fence = 0;
        end
        
        
        % Calculate the rotation matrix for a change of reference frame
        % aligned with the vector [u,w]'.  
        
%         
%         % Calculating the surface coordinates
%         x = xpri;
%         z = real(tmax*c/0.2*(a0*sqrt(x/c) + a1*(x/c) + a2*(x/c).^2 + a3*(x/c).^3 + a4*(x/c).^4));
% 
%         if z0 >= 0
%             z = z-zt(1);
%             dir = 1;
%         else
%             z = -z-zb(1);
%             dir = -1;
%         end
% 
%         f = dir*(zpri - z);
%         
%         if (z0 ~= 0) && f > 0
%             fence = 0;
%         else
%             fence = 1;
%         end
        
        if fence == 1
            
            % Calculating small region parallel to the vector [u,w]'.
            p0_plus = [x0;z0] + Lmax*n';
            p0_minus = [x0;z0] - Lmax*n';

            % Combine xp and zp.
            pp = zeros(2*(Npan + 1),1);
            pp(1:2:end-1) = xp;
            pp(2:2:end) = zp;

            % Rotating body points into coordinate system [u,w]'.
            P0 = kron(eye(Npan + 1),R0)*(pp - kron(ones(Npan + 1,1),[x0;z0]));
            P0x = P0(1:2:end-1);
            Pplus = kron(eye(Npan + 1),R0)*(pp - kron(ones(Npan + 1,1),p0_plus));
            Pminus = kron(eye(Npan + 1),R0)*(pp - kron(ones(Npan + 1,1),p0_minus));

            % Finding points that are less than p0_plus and greater than
            % p0_minus.
            IndRed = find((Pplus(2:2:end) < 0) & (Pminus(2:2:end) > 0));

            % Defining a reduced set of body points.
            xp_red = xp(IndRed);
            zp_red = zp(IndRed);
            vnp_red = vnp(IndRed,:);
            P0x_red = P0x(IndRed);

            % Defining a vector from [x0,z0] to each body point.
            r0 = [xp_red - x0 zp_red - z0];
            D0 = sqrt(r0(:,1).^2 + r0(:,2).^2);
            r0_hat = diag(1./D0)*r0;

            %  The first surface is the one with n dot r0_hat < 0.  
            Fsurf = sum(vnp_red.*r0_hat,2);

            % Reducing the set of body points further to the first surface
            % only.
            IndRed2 = find(Fsurf < 0);

            xp_red2 = xp_red(IndRed2);
            zp_red2 = zp_red(IndRed2);
            vnp_red2 = vnp_red(IndRed2,:);
            P0x_red2 = P0x_red(IndRed2);
     
            
            % Calculating the minimum and maximum values of X' in the 
            % rotated reference frame. 
            Xmin = min(P0x_red2);
            Xmax = max(P0x_red2);
             
            % Calculating the minimum and maximum values of C necessary to
            % capture the intersection point.
            Cmin = Xmin/norm([u w])/delT;
            Cmax = Xmax/norm([u w])/delT;
            
            % Calculating a range of C.
            C = linspace(Cmin,Cmax,N)';
            
            % Initializing vectors.
            D1 = zeros(N,1); 
            D2 = zeros(N,1);
            Val1 = zeros(N,1); 
            Val2 = zeros(N,1);
            min1 = zeros(N,1); 
            min2 = zeros(N,1);
            
            for j = 1:N
                xpri = x0 + C(j)*u*delT;
                zpri = z0 + C(j)*w*delT;

                % Calculating the vector from the point (x',z') to all of the
                % body panel points. Calculating the distance and the unit vector
                % r_hat.
%                 r = [kron(xp_red2,ones(N,1)) - kron(ones(length(xp_red2),1),xpri) kron(zp_red2,ones(N,1)) - kron(ones(length(zp_red2),1),zpri)];
                r = [xp_red2 - xpri zp_red2 - zpri];
                D = sqrt(r(:,1).^2 + r(:,2).^2);
                r_hat = diag(1./D)*r;

                % Finding minimum distance points (xmin1,zmin1) and
                % (xmin2,zmin2)
                [~,Ind] = sort(D);
                min1(j) = Ind(1);
                min2(j) = Ind(2);
                D1(j) = D(Ind(1));
                D2(j) = D(Ind(2));

                % Calculating (n_hat dot r_hat)
                Val1(j) = dot(vnp_red2(Ind(1),:),r_hat(Ind(1),:));
                Val2(j) = dot(vnp_red2(Ind(2),:),r_hat(Ind(2),:));
            end
            
            
%             % Finding the first case where (x',z') crosses the boundary of
%             % the body.  
%             CInd = find(Val1 > 0,1,'first');
%             
%             % Finding the nearest surface points.
%             xp1 = xp_red2(min1(CInd))
%             zp1 = zp_red2(min1(CInd))
%             
%             xp2 = xp_red2(min2(CInd))
%             zp2 = zp_red2(min2(CInd))

            % Combine xp_red2 and zp_red2.
            pp_red2 = zeros(2*length(xp_red2),1);
            pp_red2(1:2:end-1) = xp_red2;
            pp_red2(2:2:end) = zp_red2;
            
            % Rotating the nearest surface points into the frame of
            % reference of the vector [u,w]'.
            P0_red2 = kron(eye(length(xp_red2)),R0)*(pp_red2 - kron(ones(length(xp_red2),1),[x0;z0]));

            P0_X = P0_red2(1:2:end-1);
            P0_Z = P0_red2(2:2:end);
            
%             Xp1 = R0*([xp1;zp1] - [x0;z0])
%             Xp2 = R0*([xp2;zp2] - [x0;z0])
            
%             Xp = [Xp1(2);Xp2(2)]
%             Zp = -[Xp1(1);Xp2(1)]
            Xp = P0_Z;
            Zp = -P0_X;
            
            % Calculating (X*,Z*).
            Xstar = -interp1(Xp,Zp,0,'pchip');
            Zstar = 0;
            
            % Rotating points back to global reference frame.
            pstar = R0'*[Xstar; Zstar] + [x0;z0];
            xstar = pstar(1);
            zstar = pstar(2);
            
            % Find the two nearest points to the intersection point.
            r = [xp_red2 - xstar zp_red2 - zstar];
            D = sqrt(r(:,1).^2 + r(:,2).^2);
            [~,Ind] = sort(D);
            
            % Calculating the body normal at the intersection point.
            vnbod = vnp_red2(Ind(1),:)/2 + vnp_red2(Ind(2),:)/2;       
            
            % Calculating C value.
            Cvalz = (xstar - x0)/u/delT;
            
            Cval = Cvalz;
            


%             % Calculating the surface coordinates
%             x = xpri;
%             z = real(tmax*c/0.2*(a0*sqrt(x/c) + a1*(x/c) + a2*(x/c).^2 + a3*(x/c).^3 + a4*(x/c).^4));
% 
%             if z0 >= 0
%                 z = z-zt(1);
%             else
%                 z = -z-zb(1);
%             end

            

%             f = zpri - z;
% 
% 
%             % Interpolating the intersection location
%             Ci = interp1(f,C,0,'pchip');

            xstar = x0 + Cval*u*delT;
            zstar = z0 + Cval*w*delT;


%             % Calculating the body normal at the intersection point.
%             tbodvec = [xp1 - xp2 zp1 - zp2];
%             vtbod = tbodvec/norm(tbodvec);
%             vnbod = [-vtbod(2) vtbod(1)];
            
            
%             xplus = xstar + 1e-5;
%             xminus = xstar - 1e-5;
% 
%             zplus = tmax*c/0.2*(a0*sqrt(xplus/c) + a1*(xplus/c) + a2*(xplus/c).^2 + a3*(xplus/c).^3 + a4*(xplus/c).^4);
%             zminus = tmax*c/0.2*(a0*sqrt(xminus/c) + a1*(xminus/c) + a2*(xminus/c).^2 + a3*(xminus/c).^3 + a4*(xminus/c).^4);
% 
%             if zstar > 0
%                 zplus = zplus-zt(1);
%                 zminus = zminus-zt(1);
%             elseif zstar == 0
%                 zplus = zplus-zt(1);
%                 zminus = -zplus-zb(1);
%             else
%                 zplus = -zplus-zb(1);
%                 zminus = -zminus-zb(1);
%             end

%             vtbod = real([sign(z0)*(xplus - xminus); sign(z0)*(zplus - zminus)]./norm([(xplus - xminus); (zplus - zminus)]));
%             vnbod = real([-vtbod(2); vtbod(1)]);

            % Moving intersection point outward by epsilon
            xstar_n(i) = xstar + epBod*vnbod(1);
            zstar_n(i) = zstar + epBod*vnbod(2);

            % % Plotting surface, possible coordinates (x',z'), the intersection
            % % coordinate (x*, z*) and the outwardly moved intersection coordinate.
            % fighand = figure;
            % hold on
            % axis equal
            % plot(x_c*c,zt,'-k','linewidth',2)
            % plot(x_c*c,zb,'-k','linewidth',2)
            % plot(x0,z0,'ob','linewidth',2)
            % plot(xstar,zstar,'xg','linewidth',2)
            % plot(xstar_n,zstar_n,'xb','linewidth',2)
            % plot(xpri(end),zpri(end),'or','linewidth',2)
            % plot(xpri,zpri,'-b','linewidth',1)
            % 
            % xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
            % ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
        end
    end
    
    xstar_n = xstar_n + x_b;
    zstar_n = zstar_n + z_b;
end
