% clear
% clc
% 
% %% Importing Data
% load('MR_FlowVisTest.mat');


%% Calculating Off-body Flowfield

% for i_t0 = 1:Ncyc*Nstep + 1;
i_t0 = Ncyc*Nstep + 1;

    if i_t0 == 1;
        d = 1;
    else
        d = i_t0-1;
    end

    t0 = (i_t0 - 1)*delT;

%     yplane = 0.55*(1+b_b)*b;
    yplane = 3;

%     % Calculating shift in panel positions with the freestream velocity. 
%     delpc = kron(kron([1 1 1 1]',[x_b(i_t0) y_b(i_t0) z_b(i_t0)]' + D*[0 0 1]'),ones(1,2*Npanels));
%     delpcTE = kron(kron([1 1 1 1]',[x_b(i_t0) y_b(i_t0) z_b(i_t0)]' + D*[0 0 1]'),ones(1,Nspanels));
%     delcpts = kron([x_b(i_t0) y_b(i_t0) z_b(i_t0)]' + D*[0 0 1]',ones(1,2*Npanels));
%     delcptsTE = kron([x_b(i_t0) y_b(i_t0) z_b(i_t0)]' + D*[0 0 1]',ones(1,Nspanels));
%     delP(:,:,1) = x_b(i_t0)*ones(2*Nx - 1,Ny);
%     delP(:,:,2) = y_b(i_t0)*ones(2*Nx - 1,Ny);
%     delP(:,:,3) = z_b(i_t0)*ones(2*Nx - 1,Ny) + D*ones(2*Nx - 1,Ny);
%     pcTEtemp = pcTE;

    % Calculating new positions of the fin surface.    
    %     Atip = Atipvec(i_t);
%         [pc,Vc,cpts,vc,vn,vt,S,P,C1,dl,dm,dm_hat] = KinematicsHeavePitch(thick,c,Ny,Nx,yf,xf,t,f,phi,LE,alpha_max,h_c,c_r);
%       [pc,Vc,cpts,vc,vn,vt,S,P,C1,dl,dm,dm_hat] = Kinematics(thick,c,Ny,Nx,yf,yfin,xf,b,t0,f,Kx,Ky,n,phi,LE,bod,b_b);
    %     [pc,Vc,cpts,vc,vn,vt,S,P,C1,dl,dm,dm_hat] = KinematicsEllp(thick,c,Ny,Nx,yf,xf,b,t,f,Kx,Ky,phi,LE,Atip);

    % Calculating the new trailing edge panels.
%     [pcTE,cptsTE,vcTE,vnTE,vtTE,STE,dlTE,dmTE] = TEPanel(pc,Q0,delT,Ncpanels,Nspanels);    
% 
%     % Applying freestream velocity to panels.
%     pc = pc + delpc;
%     pcTE = pcTE + delpcTE;
%     cpts = cpts + delcpts;
%     cptsTE = cptsTE + delcptsTE;
%     P = P + delP;

    wvec = 1:(d-1)*Nspanels + Nspanels;


%     [MC,LE,TE,thick,c,xf] = Geometry(yplane,b,c_r,c_ms,c_t,x_ms,x_3qs,x_t,Nx,1,tmax_b,tmax_f,ybod,yfin,c_b);
%     [MC_ellp,LE_ellp,TE_ellp,thick_ellp,c_ellp,xf_ellp] = GeometryElliptical(yplane,b,c_r,Nx,tmax_f);
%     [P_ellp] = KinematicsEllpAirfoil(thick_ellp,c_ellp,1,Nx,yplane,xf_ellp,b,t0,f,Kx,Ky,phi,LE_ellp,Atip);
%     delP_ellp(:,1) = x_b(i_t0)*ones(2*Nx - 1,1);
%     delP_ellp(:,2) = y_b(i_t0)*ones(2*Nx - 1,1);
%     delP_ellp(:,3) = z_b(i_t0)*ones(2*Nx - 1,1) + D*ones(2*Nx - 1,1);
%     P_ellp = P_ellp + delP_ellp;
    P_ellp = [];

    tic
    nx = 30;
    ny = nx;
    nz = nx;
    SC = 1;
    ep = 1e-2;
    epSC = ep;
%     ep = min([Qinf*delT b/Nspanels])*5e-1;

    U = -Q0(1,i_t0);
    W = -Q0(3,i_t0);
%     U = Qinf*cos(alpha)*ones(nz,ny);
%     W = -Qinf*sin(alpha)*ones(nz,ny);
    xplane = linspace(x_b(end-1)-c_r,x_b(end-1) + 3*c_r,nx)';
%     xplane = linspace(min(P_ellp(:,1)),max(P_ellp(:,1))+ 3*c_r,nx)';
%     xf = linspace(Q0(1)*delT*(Ncyc*Nstep+1)+3/4*c_r,Q0(1)*delT*(Ncyc*Nstep+1) + 3.5*c_r,nx)';
    % yf = linspace(-4*c(1) + Q0(2)*t,c(1) + Q0(2)*t,ny)';
    % yf = linspace(0.1,0,ny)';
    % zf = linspace(c(1),-2*c(1),nz)';
    zplane = linspace(z_b(end-1) + c_r,z_b(end-1) - c_r,nz)';
    % zf = linspace(2*Atip,-2*Atip,nz)';
    % zf = 0;
    % zf = linspace(c(1) + Q0(3)*t,-c(1) + Q0(3)*t,nz)';


    % yf = linspace(-1.8*b,1.8*b,ny)';
    % yf = 0;
    % zf = 0;
    % X = dm(1)/2;
    [Xplane,Zplane] = meshgrid(xplane,zplane);
    Yplane = yplane*ones(nz,nx);
    % [Yf,Zf] = meshgrid(yf,zf);
    % [Xf,Yf] = meshgrid(xf,yf);

    X = zeros(1,nx*nz);
    Y = zeros(1,nx*nz);
    Z = zeros(1,nx*nz);

    for i = 1:nz
        vec = (i-1)*nx + 1:i*nx;
        X(1,vec) = Xplane(i,:);
        Y(1,vec) = Yplane(i,:);
        Z(1,vec) = Zplane(i,:);
    end

    ubody = zeros(nz,nx,2);
    vbody = zeros(nz,nx,2);
    wbody = zeros(nz,nx,2);
    uTE = zeros(nz,nx,2);
    vTE = zeros(nz,nx,2);
    wTE = zeros(nz,nx,2);
    uw = zeros(nz,nx,2);
    vw = zeros(nz,nx,2);
    ww = zeros(nz,nx,2);
    ud = zeros(nz,nx);
    vd = zeros(nz,nx);
    wd = zeros(nz,nx);
    us = zeros(nz,nx);
    vs = zeros(nz,nx);
    ws = zeros(nz,nx);
    udI = zeros(nz,nx);
    vdI = zeros(nz,nx);
    wdI = zeros(nz,nx);
    usI = zeros(nz,nx);
    vsI = zeros(nz,nx);
    wsI = zeros(nz,nx);


    % Calculating induced velocity from body panels.
    [udt,vdt,wdt,ust,vst,wst] = InfluenceVMatrix(X,Y,Z,pc,cpts,vc,vn,vt,mu(:,i_t0),sig,dl,dm,S,ep,epSC,0);
    [udIt,vdIt,wdIt,usIt,vsIt,wsIt] = InfluenceVMatrix(X,Y,Z,pcI,cptsI,vcI,vnI,vtI,mu(:,i_t0),sig,dlI,dmI,S,ep,epSC,0);

    for i = 1:nz
        vec = (i-1)*nx + 1:i*nx;
        ud(i,:) = udt(1,vec);
        vd(i,:) = vdt(1,vec);
        wd(i,:) = wdt(1,vec);
        us(i,:) = ust(1,vec);
        vs(i,:) = vst(1,vec);
        ws(i,:) = wst(1,vec);

        udI(i,:) = udIt(1,vec);
        vdI(i,:) = vdIt(1,vec);
        wdI(i,:) = wdIt(1,vec);
        usI(i,:) = usIt(1,vec);
        vsI(i,:) = vsIt(1,vec);
        wsI(i,:) = wsIt(1,vec);
    end

    ud2 = zeros(nz,nx);
    vd2 = zeros(nz,nx);
    wd2 = zeros(nz,nx);
    us2 = zeros(nz,nx);
    vs2 = zeros(nz,nx);
    ws2 = zeros(nz,nx);
    udI2 = zeros(nz,nx);
    vdI2 = zeros(nz,nx);
    wdI2 = zeros(nz,nx);
    usI2 = zeros(nz,nx);
    vsI2 = zeros(nz,nx);
    wsI2 = zeros(nz,nx);

    if grd == 1
        [ud2t,vd2t,wd2t,us2t,vs2t,ws2t] = InfluenceVMatrix(X,Y,Z,pc2,cpts2,vc2,vn2,vt2,mu(:,i_t0),sig,dl2,dm2,S2,ep,epSC,0);
        [udI2t,vdI2t,wdI2t,usI2t,vsI2t,wsI2t] = InfluenceVMatrix(X,Y,Z,pcI2,cptsI2,vcI2,vnI2,vtI2,mu(:,i_t0),sig,dlI2,dmI2,SI2,ep,epSC,0);

        for i = 1:nz
            vec = (i-1)*nx + 1:i*nx;
            ud2(i,:) = ud2t(1,vec);
            vd2(i,:) = vd2t(1,vec);
            wd2(i,:) = wd2t(1,vec);
            us2(i,:) = us2t(1,vec);
            vs2(i,:) = vs2t(1,vec);
            ws2(i,:) = ws2t(1,vec);

            udI2(i,:) = udI2t(1,vec);
            vdI2(i,:) = vdI2t(1,vec);
            wdI2(i,:) = wdI2t(1,vec);
            usI2(i,:) = usI2t(1,vec);
            vsI2(i,:) = vsI2t(1,vec);
            wsI2(i,:) = wsI2t(1,vec);
        end

    end

    ubody(:,:,1) = ud + udI + ud2 + udI2; 
    ubody(:,:,2) = us + usI + us2 + usI2;
    vbody(:,:,1) = vd + vdI + vd2 + vdI2;
    vbody(:,:,2) = vs + vsI + vs2 + vsI2;
    wbody(:,:,1) = wd + wdI + wd2 + wdI2;
    wbody(:,:,2) = ws + wsI + ws2 + wsI2;

    % Calculating induced velocity from trailing edge panels.
    [udt,vdt,wdt,ust,vst,wst] = InfluenceVMatrix(X,Y,Z,pcTE,cptsTE,vcTE,vnTE,vtTE,muTE(:,i_t0),0*muTE(:,i_t0),dlTE,dmTE,STE,ep,epSC,0);
    [udIt,vdIt,wdIt,usIt,vsIt,wsIt] = InfluenceVMatrix(X,Y,Z,pcTEI,cptsTEI,vcTEI,vnTEI,vtTEI,muTE(:,i_t0),0*muTE(:,i_t0),dlTEI,dmTEI,STEI,ep,epSC,0);

    for i = 1:nz
        vec = (i-1)*nx + 1:i*nx;
        ud(i,:) = udt(1,vec);
        vd(i,:) = vdt(1,vec);
        wd(i,:) = wdt(1,vec);
        us(i,:) = ust(1,vec);
        vs(i,:) = vst(1,vec);
        ws(i,:) = wst(1,vec);

        udI(i,:) = udIt(1,vec);
        vdI(i,:) = vdIt(1,vec);
        wdI(i,:) = wdIt(1,vec);
        usI(i,:) = usIt(1,vec);
        vsI(i,:) = vsIt(1,vec);
        wsI(i,:) = wsIt(1,vec);
    end

    ud2 = zeros(nz,nx);
    vd2 = zeros(nz,nx);
    wd2 = zeros(nz,nx);
    us2 = zeros(nz,nx);
    vs2 = zeros(nz,nx);
    ws2 = zeros(nz,nx);
    udI2 = zeros(nz,nx);
    vdI2 = zeros(nz,nx);
    wdI2 = zeros(nz,nx);
    usI2 = zeros(nz,nx);
    vsI2 = zeros(nz,nx);
    wsI2 = zeros(nz,nx);

    if grd == 1
        [ud2t,vd2t,wd2t,us2t,vs2t,ws2t] = InfluenceVMatrix(X,Y,Z,pcTE2,cptsTE2,vcTE2,vnTE2,vtTE2,muTE(:,i_t0),0*muTE(:,i_t0),dlTE2,dmTE2,STE2,ep,epSC,0);
        [udI2t,vdI2t,wdI2t,usI2t,vsI2t,wsI2t] = InfluenceVMatrix(X,Y,Z,pcTEI2,cptsTEI2,vcTEI2,vnTEI2,vtTEI2,muTE(:,i_t0),0*muTE(:,i_t0),dlTEI2,dmTEI2,STEI2,ep,epSC,0);

        for i = 1:nz
            vec = (i-1)*nx + 1:i*nx;
            ud2(i,:) = ud2t(1,vec);
            vd2(i,:) = vd2t(1,vec);
            wd2(i,:) = wd2t(1,vec);
            us2(i,:) = us2t(1,vec);
            vs2(i,:) = vs2t(1,vec);
            ws2(i,:) = ws2t(1,vec);

            udI2(i,:) = udI2t(1,vec);
            vdI2(i,:) = vdI2t(1,vec);
            wdI2(i,:) = wdI2t(1,vec);
            usI2(i,:) = usI2t(1,vec);
            vsI2(i,:) = vsI2t(1,vec);
            wsI2(i,:) = wsI2t(1,vec);
        end
    end

    uTE(:,:,1) = ud + udI + ud2 + udI2;  
    uTE(:,:,2) = us + usI + us2 + usI2;
    vTE(:,:,1) = vd + vdI + vd2 + vdI2;
    vTE(:,:,2) = vs + vsI + vs2 + vsI2;
    wTE(:,:,1) = wd + wdI + wd2 + wdI2;
    wTE(:,:,2) = ws + wsI + ws2 + wsI2;

    % Calculating induced velocity from wake panels.  
    [udt,vdt,wdt,ust,vst,wst] = InfluenceVMatrix(X,Y,Z,pcw(:,wvec),cptsw(:,wvec),vcw(:,wvec),vnw(:,wvec),vtw(:,wvec),muW(wvec),0*muW(wvec),dlw(wvec),dmw(wvec),Sw(wvec),ep,epSC,SC);
    [udIt,vdIt,wdIt,usIt,vsIt,wsIt] = InfluenceVMatrix(X,Y,Z,pcwI(:,wvec),cptswI(:,wvec),vcwI(:,wvec),vnwI(:,wvec),vtwI(:,wvec),muW(wvec),0*muW(wvec),dlwI(wvec),dmwI(wvec),SwI(wvec),ep,epSC,SC);

    for i = 1:nz
        vec = (i-1)*nx + 1:i*nx;
        ud(i,:) = udt(1,vec);
        vd(i,:) = vdt(1,vec);
        wd(i,:) = wdt(1,vec);
        us(i,:) = ust(1,vec);
        vs(i,:) = vst(1,vec);
        ws(i,:) = wst(1,vec);

        udI(i,:) = udIt(1,vec);
        vdI(i,:) = vdIt(1,vec);
        wdI(i,:) = wdIt(1,vec);
        usI(i,:) = usIt(1,vec);
        vsI(i,:) = vsIt(1,vec);
        wsI(i,:) = wsIt(1,vec);
    end

    ud2 = zeros(nz,nx);
    vd2 = zeros(nz,nx);
    wd2 = zeros(nz,nx);
    us2 = zeros(nz,nx);
    vs2 = zeros(nz,nx);
    ws2 = zeros(nz,nx);
    udI2 = zeros(nz,nx);
    vdI2 = zeros(nz,nx);
    wdI2 = zeros(nz,nx);
    usI2 = zeros(nz,nx);
    vsI2 = zeros(nz,nx);
    wsI2 = zeros(nz,nx);

    if grd == 1
        [ud2t,vd2t,wd2t,us2t,vs2t,ws2t] = InfluenceVMatrix(X,Y,Z,pcw2(:,wvec),cptsw2(:,wvec),vcw2(:,wvec),vnw2(:,wvec),vtw2(:,wvec),muW(wvec),0*muW(wvec),dlw2(wvec),dmw2(wvec),Sw2(wvec),ep,epSC,SC);
        [udI2t,vdI2t,wdI2t,usI2t,vsI2t,wsI2t] = InfluenceVMatrix(X,Y,Z,pcwI2(:,wvec),cptswI2(:,wvec),vcwI2(:,wvec),vnwI2(:,wvec),vtwI2(:,wvec),muW(wvec),0*muW(wvec),dlwI2(wvec),dmwI2(wvec),SwI2(wvec),ep,epSC,SC);

        for i = 1:nz
            vec = (i-1)*nx + 1:i*nx;
            ud2(i,:) = ud2t(1,vec);
            vd2(i,:) = vd2t(1,vec);
            wd2(i,:) = wd2t(1,vec);
            us2(i,:) = us2t(1,vec);
            vs2(i,:) = vs2t(1,vec);
            ws2(i,:) = ws2t(1,vec);

            udI2(i,:) = udI2t(1,vec);
            vdI2(i,:) = vdI2t(1,vec);
            wdI2(i,:) = wdI2t(1,vec);
            usI2(i,:) = usI2t(1,vec);
            vsI2(i,:) = vsI2t(1,vec);
            wsI2(i,:) = wsI2t(1,vec);
        end
    end

    uw(:,:,1) = ud + udI + ud2 + udI2; 
    uw(:,:,2) = us + usI + us2 + usI2;
    vw(:,:,1) = vd + vdI + vd2 + vdI2;
    vw(:,:,2) = vs + vsI + vs2 + vsI2;
    ww(:,:,1) = wd + wdI + wd2 + wdI2;
    ww(:,:,2) = ws + wsI + ws2 + wsI2;

    % Calculating total induced velocity from doublet and source panels.

    u_d = ubody(:,:,1) + uTE(:,:,1) + uw(:,:,1);
    v_d = vbody(:,:,1) + vTE(:,:,1) + vw(:,:,1);
    w_d = wbody(:,:,1) + wTE(:,:,1) + ww(:,:,1);

    u_s = ubody(:,:,2) + uTE(:,:,2) + uw(:,:,2);
    v_s = vbody(:,:,2) + vTE(:,:,2) + vw(:,:,2);
    w_s = wbody(:,:,2) + wTE(:,:,2) + ww(:,:,2);

    u = u_d + u_s;
    v = v_d + v_s;
    w = w_d + w_s;

    U = U + u;
    W = W + w;

    % Calculating vorticity:  
    % omega_y = du/dz - dw/dx, omega_z = dv/dx - du/dy.
    omega_y = (u(3:end,2:end-1) - u(1:end-2,2:end-1))./(Zplane(3:end,2:end-1) - Zplane(1:end-2,2:end-1)) - (w(2:end-1,3:end) - w(2:end-1,1:end-2))./(Xplane(2:end-1,3:end) - Xplane(2:end-1,1:end-2));
    % omega_z = (v(2:end-1,3:end) - v(2:end-1,1:end-2))./(Xplane(2:end-1,3:end) - Xplane(2:end-1,1:end-2)) - (u(3:end,2:end-1) - u(1:end-2,2:end-1))./(Yplane(3:end,2:end-1) - Yplane(1:end-2,2:end-1));
    Xstar = (Xplane(2:end-1,3:end) + Xplane(2:end-1,1:end-2))/2;
    % Ystar = (Yplane(3:end,2:end-1) + Yplane(1:end-2,2:end-1))/2;
    Zstar = (Zplane(3:end,2:end-1) + Zplane(1:end-2,2:end-1))/2;

    [rows,cols] = size(Xplane);
    toc


    %% Plotting flowfields (X-Z plane).

%     f1 = figure
%     hold on
%     plot3(P(:,:,1) - Q0(1)*t0,P(:,:,2),P(:,:,3),'.k','linewidth',2,'markersize',12);
%     % plot3(P_ellp(:,1),P_ellp(:,2),P_ellp(:,3),'.-k','linewidth',2,'markersize',12);
%     % ptch = patch(P_ellp(:,1),P_ellp(:,2),P_ellp(:,3),'k');
%     % plot((P(:,1,2))/c(1),(P(:,1,3))/c(1),'-k','linewidth',1.5);
%     % streamline(-(Yf - Q0(2)*t)/c(1),(Zf - Q0(3)*t)/c(1),-V,W,-(Yf(:,end) - Q0(2)*t)/c(1),(Zf(:,end) - Q0(3)*t)/c(1));
%     % quiver(Yf,Zf,-v,w)
%     % pcolor(Xstar,yf*ones(rows,cols),Zstar,omega_y)
%     quiver3((Xplane - Q0(1)*t0),yplane*ones(rows,cols),Zplane,u,0*w,w,'k')
%     pcolor(Xstar,Zstar,omega_y)
%     shading interp
%     view(0,0)
%     axis equal
%     % axis([-0.2 0.4 0 b 0 0.2])
%     xlabel('$$x/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%     ylabel('$$z/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    
    
    f2 = figure;
    hold on
%     pcolor((Xstar - max(P_ellp(:,1)))/c_r,(Zstar)/c_r,-omega_y)
%     shading interp
%     plot((P_ellp(:,1) - max(P_ellp(:,1)))/c_r,P_ellp(:,3)/c_r,'.-k','linewidth',2,'markersize',12);
    % ptch = patch(P_ellp(:,1),P_ellp(:,3),'k');
%     quiver(Xplane,Zplane,U,W,'b')
    num = 101
    startx = Xplane(1,1)*ones(1,num);
    startz = linspace(Zplane(1,1),Zplane(end,1)-3*c_r,num);
%     streamline(-(Yf - Q0(2)*t)/c(1),(Zf - Q0(3)*t)/c(1),-V,W,-(Yf(:,end) - Q0(2)*t)/c(1),(Zf(:,end) - Q0(3)*t)/c(1));
    streamline(Xplane,Zplane,U,W,startx,startz) 
%     colorbar
%     caxis(1*[-1 1])
    axis equal
    % axis([-5 -1 -b b -b b])
    xlabel('$$x/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    ylabel('$$z/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

    [xp,Cp] = PanelCode2D_Cp(-alpha,100)
    
    % figure
    % hold on
    % pcolor(Xstar,Zstar,-omega_y)
    % shading interp
    % plot(P_ellp(:,1),P_ellp(:,3),'.-k','linewidth',2,'markersize',12);
    % % ptch = patch(P_ellp(:,1),P_ellp(:,3),'k');
    % % quiver(Xf,Zf,u_d + u_s - Vt(1,end),w_s + w_d,'k')
    % colorbar
    % % caxis([-2 2])
    % axis equal
    % % axis([-5 -1 -b b -b b])
    % xlabel('$$x/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    % ylabel('$$z/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

    % figure
    % hold on
    % quiver((Xf - max(P_ellp(:,1)))/c_r,Zf/c_r,u_d,w_d,'k')
    % axis equal
    % xlabel('$$x/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    % ylabel('$$z/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

    % figure
    % hold on
    % quiver((Xf - max(P_ellp(:,1)))/c_r,Zf/c_r,u_s,w_s,'k')
    % axis equal
    % xlabel('$$x/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    % ylabel('$$z/c_r$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

    % %% Plotting flowfields (X-Y plane).
    % figure
    % hold on
    % pcolor(Xstar,Ystar,omega_z)
    % shading interp
    % 
    % % plot3(P(:,:,1),P(:,:,2),P(:,:,3),'.k','markersize',10);
    % % plot3(P(:,:,1),-P(:,:,2),P(:,:,3),'.k','markersize',10);
    % PlotPanels(t,pc,pcI,Vc,pcTE,pcTEI,pcw,pcwI,C1,cpts,cptsI,cptsTE,cptsTEI,cptsw,cptswI,vc,vcI,vn,vnI,vt,vtI,vcTE,vcTEI,vnTE,vnTEI,vtTE,vtTEI,vcw,vcwI,vnw,vnwI,vtw,vtwI,Npanels,Nspanels,Cp(i_t(end),:))
    % if grd == 1
    %     PlotPanels(t,pc2,pcI2,Vc,pcTE2,pcTEI2,pcw2,pcwI2,C1,cpts2,cptsI2,cptsTE2,cptsTEI2,cptsw2,cptswI2,vc2,vcI2,vn2,vnI2,vt2,vtI2,vcTE2,vcTEI2,vnTE2,vnTEI2,vtTE2,vtTEI2,vcw2,vcwI2,vnw2,vnwI2,vtw2,vtwI2,Npanels,Nspanels,Cp(i_t,:))
    %     PlotGrd(Ncyc,Nstep,delT,Q0,c_b,b,c_r,C1);
    % end
    % quiver3(Xf,Yf,zf*ones(rows,cols),u,v,0*v,'k')
    % colorbar
    % caxis([-1 1])
    % view(2)
    % axis equal
    % % axis([-12.5 -7.5 -1.8*b 1.8*b Ncyc*Nstep*delT*Q0(3)-1.5*b 1.5*b])
    % % axis([(Ncyc*Nstep + 2)*delT*Q0(1) - (c_b*b - c_r)/2*1.05 1.05*c_b*b -1.6*b 1.6*b Ncyc*Nstep*delT*Q0(3)-1.5*b 1.5*b])
    % xlabel('$$x, m$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    % ylabel('$$y, m$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

    % figure
    % hold on
    % PlotPanels(t,pc,pcI,Vc,pcTE,pcTEI,pcw,pcwI,C1,cpts,cptsI,cptsTE,cptsTEI,cptsw,cptswI,vc,vcI,vn,vnI,vt,vtI,vcTE,vcTEI,vnTE,vnTEI,vtTE,vtTEI,vcw,vcwI,vnw,vnwI,vtw,vtwI,Npanels,Nspanels,Cp(i_t(end),:))
    % if grd == 1
    %     PlotPanels(t,pc2,pcI2,Vc,pcTE2,pcTEI2,pcw2,pcwI2,C1,cpts2,cptsI2,cptsTE2,cptsTEI2,cptsw2,cptswI2,vc2,vcI2,vn2,vnI2,vt2,vtI2,vcTE2,vcTEI2,vnTE2,vnTEI2,vtTE2,vtTEI2,vcw2,vcwI2,vnw2,vnwI2,vtw2,vtwI2,Npanels,Nspanels,Cp(i_t,:))
    %     PlotGrd(Ncyc,Nstep,delT,Q0,c_b,b,c_r,C1);
    % end
    % quiver3(Xf,Yf,zf*ones(rows,cols),u,v,0*v,'k')
    % view(2)
    % axis equal
    % axis([-6 -1 -1.55*b 1.55*b -b b])
    % xlabel('$$x, m$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    % ylabel('$$y, m$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    % 
    % figure
    % hold on
    % PlotPanels(t,pc,pcI,Vc,pcTE,pcTEI,pcw,pcwI,C1,cpts,cptsI,cptsTE,cptsTEI,cptsw,cptswI,vc,vcI,vn,vnI,vt,vtI,vcTE,vcTEI,vnTE,vnTEI,vtTE,vtTEI,vcw,vcwI,vnw,vnwI,vtw,vtwI,Npanels,Nspanels,Cp(i_t(end),:))
    % if grd == 1
    %     PlotPanels(t,pc2,pcI2,Vc,pcTE2,pcTEI2,pcw2,pcwI2,C1,cpts2,cptsI2,cptsTE2,cptsTEI2,cptsw2,cptswI2,vc2,vcI2,vn2,vnI2,vt2,vtI2,vcTE2,vcTEI2,vnTE2,vnTEI2,vtTE2,vtTEI2,vcw2,vcwI2,vnw2,vnwI2,vtw2,vtwI2,Npanels,Nspanels,Cp(i_t,:))
    %     PlotGrd(Ncyc,Nstep,delT,Q0,c_b,b,c_r,C1);
    % end
    % quiver3(Xf,Yf,zf*ones(rows,cols),u,v,0*v,'k')
    % view(2)
    % axis equal
    % axis([-6 -1 -1.55*b 1.55*b -b b])
    % xlabel('$$x, m$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    % ylabel('$$y, m$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')


%     Flow(i_t0) = getframe;
%     close(f2)
% end


% movie(Flow,5,3/delT)