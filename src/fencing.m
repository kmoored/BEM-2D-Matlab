
function [xstar_n,zstar_n,fence] = fencing(c,Vp,xpts,zpts,up,wp,xp,zp,vn,x_b,z_b,delT,epB)

    N = 50;
    
    xpts = xpts - x_b;
    zpts = zpts - z_b;
    
    xp = xp - x_b;
    zp = zp - z_b;
    
    % Defining the top and bottom surfaces
    Npan = length(xp) - 1;
    
    % Defining the outward normal vectors at each panel endpoint.
    vnp = zeros(length(xp),2);
    vnp(1,:) = vn(1,:)/2 + vn(end,:)/2;
    vnp(end,:) = vn(1,:)/2 + vn(end,:)/2;
    vnp(2:end-1,:) = vn(1:end-1,:)/2 + vn(2:end,:)/2;
    
    xp = xp + epB*vnp(:,1);
    zp = zp + epB*vnp(:,2);
    
    L = sqrt((xp(2:end) - xp(1:end-1)).^2 + (zp(2:end) - zp(1:end-1)).^2);
    Lmax = max(L);
    
    xstar_n = xpts;
    zstar_n = zpts;
    
    xpri_all = xpts + up*delT;
    zpri_all = zpts + wp*delT;
    
    
%     % Plotting surface, possible coordinates (x',z'), the intersection
%     % coordinate (x*, z*) and the outwardly moved intersection coordinate.
%     fighand = figure;
%     hold on
%     axis equal
%     plot(xp,zp,'-k','linewidth',2)
%     plot(xpts,zpts,'ob','linewidth',2)
%     plot(xpri_all,zpri_all,'xg','linewidth',2)
%     for i = 1:length(xpts)
%         plot([xpts(i);xpri_all(i)],[zpts(i) zpri_all(i)],'-b','linewidth',1)
%     end
% 
%     xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%     ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

    
    fence = zeros(1,length(xpts));
    fence2 = zeros(1,length(xpts));

    for i = 1:length(xpts)
        
        % Extracting panel starting coordinates
        x0 = xpts(i);
        z0 = zpts(i);
        
        % Extracting velocity field at starting coordinates
        u = up(i);
        w = wp(i);
        
        % Finding xinf anf zinf 10 chords away from the body along the
        % vector [u,w]'.
        vec = [u w]/norm([u w]);
        xinf = x0 - 10*c*vec(1);
        zinf = z0 - 10*c*vec(2); 
        
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
        r_hat = [D.\r(:,1) D.\r(:,2)];
        
        % Finding minimum distance point (xmin,zmin)
        [~,indmin] = min(D);
        
        % Calculating (n_hat dot r_hat)
        Val = sum(vnp(indmin,:).*r_hat(indmin,:));
                
        % Calculating the vector from the point (x0,z0) to all of the
        % body panel points. Calculating the distance and the unit vector
        % r0_hat.
        r0 = [xp - x0 zp - z0];
        D0 = sqrt(r0(:,1).^2 + r0(:,2).^2);
        r0_hat = [D0.\r0(:,1) D0.\r0(:,2)];

        
        % Finding minimum distance point (xmin,zmin)
        [~,indmin0] = min(D0);
        
        % Calculating (n_hat dot r_hat)
        ValIn = sum(vnp(indmin0,:).*r0_hat(indmin0,:));
                
        
        
        
        
        
        % Determining whether the point (x',z') is inside or outside of the
        % body.
        if Val > 0
            fence(i) = 1;
        else
            fence(i) = 0;
        end
        
        if ValIn > 0 && u == 0 && w == 0
            fence2(i) = 1;
        else
            fence2(i) = 0;
        end
        
        
        
        
        
        
        
        
        
        % Applying a fencing algorithm for points that start in the body.
        if fence2(i) == 1
            
            % Defining the nearest surface point neighbors.
            if indmin0 == 1 
                xs0 = [xp(Npan) xp(1) xp(2)]';
                zs0 = [zp(Npan) zp(1) zp(2)]';
            elseif indmin == Npan + 1
                xs0 = [xp(Npan) xp(1) xp(2)]';
                zs0 = [zp(Npan) zp(1) zp(2)]';
            else
                xs0 = [xp(indmin0 - 1) xp(indmin0) xp(indmin0 + 1)]';        
                zs0 = [zp(indmin0 - 1) zp(indmin0) zp(indmin0 + 1)]';
            end
            
            % Calculating the vector from the point (x0,z0) to the nearest 
            % surface points. Calculating the distance and the unit vector
            % r0_hat.
            rs0 = [x0 - xs0 z0 - zs0];
            Ds0 = sqrt(rs0(:,1).^2 + rs0(:,2).^2);
            
            % Determining which neighbor is the nearest to the point
            % (x0,z0).
            [~,Ind] = sort(Ds0);
            Smin1 = Ind(1);
            Smin2 = Ind(2);
            
            % Calulating the vector from the nearest surface point to
            % (x0,z0) and from the nearest surface point to the second
            % nearest surface point.
            drs0 = rs0(Smin1,:);
            dps0 = [xs0(Smin2) zs0(Smin2)] - drs0;
            
            % Calculating the velocities of the body at the nearest surface 
            % points.
            V1 = Vp(indmin0,:);
            if (indmin0 - 2 + Smin2) == 0
                V2 = Vp(Npan + 1,:);
            elseif(indmin0 - 2 + Smin2) == Npan + 2
                V2 = Vp(1,:);
            else
                V2 = Vp(indmin0 - 2 + Smin2,:);
            end
            W1 = 1 - norm(dot(drs0,dps0))/norm(dps0);
            W2 = norm(dot(drs0,dps0))/norm(dps0);
            
            % Calculating the weighted average body velocity
            V0 = -W1*V1 - W2*V2;
            
            % Moving the internal point by -V0*delT
            xpri = x0;
            zpri = z0;
            
            x0 = x0 - V0(1)*delT;
            z0 = z0 - V0(2)*delT;
            
            
%             
%                 % Plotting surface, possible coordinates (x',z'), the intersection
%                 % coordinate (x*, z*) and the outwardly moved intersection coordinate.
%                 fighand = figure;
%                 hold on
%                 axis equal
%                 plot(xp,zp,'-k','linewidth',2)
%                 plot(xpri,zpri,'ob','linewidth',2)
%                 plot(x0,z0,'xg','linewidth',2)
%                 plot([xpri;x0],[zpri;z0],'-b','linewidth',1)
%                 
%             
%                 xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%                 ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

            
            % Beginning normal fencing algorithm.
          
            % Extracting velocity field at starting coordinates
            u = V0(1);
            w = V0(2);

            % Finding xinf anf zinf 10 chords away from the body along the
            % vector [u,w]'.
            vec = [u w]/norm([u w]);
            xinf = x0 - 10*c*vec(1);
            zinf = z0 - 10*c*vec(2);

            % Calculating the advection of the panel points (x0,z0)-->(x',z')
            

            % Calculating the tangential and normal vectors to the direction
            % [u,w]'.
            tvec = [xpri - x0 zpri - z0];
            t = tvec/norm(tvec);
            n = [-t(2) t(1)];
            R0 = [t; n];

           


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

            % Eliminating double point if point 1 or point Npan + 1 is
            % found.
            IndRed(find(IndRed == Npan + 1)) = [];
            
            % Defining a reduced set of body points.
            xp_red = xp(IndRed);
            zp_red = zp(IndRed);
            vnp_red = vnp(IndRed,:);
            P0x_red = P0x(IndRed);

            % Defining a vector from [xinf,zinf] to each body point.
            rinf = [xp_red - xinf zp_red - zinf];
            Dinf = sqrt(rinf(:,1).^2 + rinf(:,2).^2);
            rinf_hat = diag(1./Dinf)*rinf;

            %  The first surface is the one with n dot r0_hat < 0.  
            Fsurf = sum(vnp_red.*rinf_hat,2);

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
            

            % Combine xp_red2 and zp_red2.
            pp_red2 = zeros(2*length(xp_red2),1);
            pp_red2(1:2:end-1) = xp_red2;
            pp_red2(2:2:end) = zp_red2;
            
            % Rotating the nearest surface points into the frame of
            % reference of the vector [u,w]'.
            P0_red2 = kron(eye(length(xp_red2)),R0)*(pp_red2 - kron(ones(length(xp_red2),1),[x0;z0]));

            P0_X = P0_red2(1:2:end-1);
            P0_Z = P0_red2(2:2:end);
            
            Xp = P0_Z;
            Zp = -P0_X;
            
            % Calculating (X*,Z*).
            Xstar = -interp1(Xp,Zp,0,'pchip');
            Zstar = 0;
            
            % Rotating points back to the global reference frame.
            pstar = R0'*[Xstar; Zstar] + [x0;z0];
            xstar = pstar(1);
            zstar = pstar(2);
            
            % Find the two nearest points to the intersection point.
            r = [xp_red2 - xstar zp_red2 - zstar];
            D = sqrt(r(:,1).^2 + r(:,2).^2);
            [~,Ind] = sort(D);
            
            % Calculating the body normal and tangent vectors at the 
            % intersection point.
            vnbod = vnp_red2(Ind(1),:)/2 + vnp_red2(Ind(2),:)/2;  
            vtbod = [vnbod(1,2) -vnbod(1,1)];
            
             
            % Calculating C value.
            Cvalx = (xstar - x0)/u/delT;
            Cvalz = (zstar - z0)/w/delT;
            
            Cval = Cvalz;
            
            if isnan(Cvalz)
                Cval = Cvalx;
            end


            xstar = x0 + Cval*u*delT;
            zstar = z0 + Cval*w*delT;
            
            Val = dot((1 - Cval)*[u w],vtbod);


            % Moving intersection point outward by epsilon
            xstar_n(i) = xstar + epB*vnbod(1);
            zstar_n(i) = zstar + epB*vnbod(2);
            
            % Moving point tangential to the surface by the projection of
            % [u w] and vtbod.
            xstar_n(i) = xstar_n(i) + (Val*delT + 1*epB)*vtbod(1);
            zstar_n(i) = zstar_n(i) + (Val*delT + 1*epB)*vtbod(2);

%             % Plotting surface, possible coordinates (x',z'), the intersection
%             % coordinate (x*, z*) and the outwardly moved intersection coordinate.
%             fighand = figure;
%             hold on
%             axis equal
%             plot(xp,zp,'-k','linewidth',2)
%             plot(x0,z0,'ob','linewidth',2)
%             plot(xstar,zstar,'xg','linewidth',2)
%             plot(xstar_n,zstar_n,'xb','linewidth',2)
%             plot(xpri(end),zpri(end),'or','linewidth',2)
%             plot(xpri,zpri,'-b','linewidth',1)
%             
%             xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%             ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%         
%         
        
        end
        
        
        
        
        
        
        
        
        
       
        
        if fence(i) == 1 && fence2(i) ~= 1
            
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
            
            % Eliminating double point if point 1 or point Npan + 1 is
            % found.
            IndRed(find(IndRed == Npan + 1)) = [];

            % Defining a reduced set of body points.
            xp_red = xp(IndRed);
            zp_red = zp(IndRed);
            vnp_red = vnp(IndRed,:);
            P0x_red = P0x(IndRed);

            % Defining a vector from [xinf,zinf] to each body point.
            rinf = [xp_red - xinf zp_red - zinf];
            Dinf = sqrt(rinf(:,1).^2 + rinf(:,2).^2);
            rinf_hat = diag(1./Dinf)*rinf;

            %  The first surface is the one with n dot r0_hat < 0.  
            Fsurf = sum(vnp_red.*rinf_hat,2);

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
            

            % Combine xp_red2 and zp_red2.
            pp_red2 = zeros(2*length(xp_red2),1);
            pp_red2(1:2:end-1) = xp_red2;
            pp_red2(2:2:end) = zp_red2;
            
            % Rotating the nearest surface points into the frame of
            % reference of the vector [u,w]'.
            P0_red2 = kron(eye(length(xp_red2)),R0)*(pp_red2 - kron(ones(length(xp_red2),1),[x0;z0]));

            P0_X = P0_red2(1:2:end-1);
            P0_Z = P0_red2(2:2:end);
            
            Xp = P0_Z;
            Zp = -P0_X;
            
            % Calculating (X*,Z*).
            Xstar = -interp1(Xp,Zp,0,'pchip');
            Zstar = 0;
            
            % Rotating points back to the global reference frame.
            pstar = R0'*[Xstar; Zstar] + [x0;z0];
            xstar = pstar(1);
            zstar = pstar(2);
            
            % Find the two nearest points to the intersection point.
            r = [xp_red2 - xstar zp_red2 - zstar];
            D = sqrt(r(:,1).^2 + r(:,2).^2);
            [~,Ind] = sort(D);
            
            % Calculating the body normal and tangent vectors at the 
            % intersection point.
            vnbod = vnp_red2(Ind(1),:)/2 + vnp_red2(Ind(2),:)/2;  
            vtbod = [vnbod(1,2) -vnbod(1,1)];
            
            
            % Calculating C value.
            Cvalx = (xstar - x0)/u/delT;
            Cvalz = (zstar - z0)/w/delT;
            
            Cval = Cvalz;
            
            if isnan(Cvalz)
                Cval = Cvalx;
            end
           

            xstar = x0 + Cval*u*delT;
            zstar = z0 + Cval*w*delT;
            
            Val = dot((1 - Cval)*[u w],vtbod);


            % Moving intersection point outward by epsilon
            xstar_n(i) = xstar + epB*vnbod(1);
            zstar_n(i) = zstar + epB*vnbod(2);
            
            % Moving point tangential to the surface by the projection of
            % [u w] and vtbod.
            xstar_n(i) = xstar_n(i) + (Val*delT + 1*epB)*vtbod(1);
            zstar_n(i) = zstar_n(i) + (Val*delT + 1*epB)*vtbod(2);

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
            % ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname'
            % ,'TimesNewRoman')
        
            
        end
        
        
    end
    
    
    xstar_n = xstar_n + x_b;
    zstar_n = zstar_n + z_b;
end
