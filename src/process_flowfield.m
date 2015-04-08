function process_flowfield(folder,savefilename)


%% Load Processed Data
load([folder,'/Processed',savefilename,'.mat']);



%% Flowfield calculation
Flowfield = 1;
figInd = 1;
epSC = 5e-2;
SC = 1;

% Calculating the flowfield around the wing.
nx = 101;
nz = 101;
Ut = zeros(nz,nx,Nstep);
Wt = zeros(nz,nx,Nstep);

if Flowfield == 1
    for i_t = (Ncyc - 1)*Nstep+1:Ncyc*Nstep+1
%     for i_t = (Ncyc - 1)*Nstep+1:Ncyc*Nstep+1
        U = Uinf*ones(nz,nx);
        W = Winf*ones(nz,nx);

        % Flow field for ground effect calculations
        xf = linspace(-c/4 + x_b(i_t),6*c + x_b(i_t),nx)';
        zf = linspace(-3*c/2,3*c/2,nz)';

        [Xf,Zf] = meshgrid(xf,zf);

        X = zeros(1,nx*nz);
        Z = zeros(1,nx*nz);

        for i = 1:nz
            vec = (i-1)*nx + 1:i*nx;

            X(1,vec) = Xf(i,:);
            Z(1,vec) = Zf(i,:);
        end

        % Calculating flow field

        % Body contribution
        [u_st,w_st,u_bdt,w_bdt] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp(:,i_t)',zp(:,i_t)',[vt(:,i_t,1) vt(:,i_t,2)]',[vn(:,i_t,1) vn(:,i_t,2)]',epSC,SC);

        u_st_2 = 0*u_st;
        w_st_2 = 0*w_st;
        u_bdt_2 = 0*u_bdt;
        w_bdt_2 = 0*w_bdt;

        if grd == 1
            [u_st_2,w_st_2,u_bdt_2,w_bdt_2] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp(:,i_t)',-zp(:,i_t)',[vt(:,i_t,1) -vt(:,i_t,2)]',[vn(:,i_t,1) -vn(:,i_t,2)]',epSC,SC);
        end

        u_s = zeros(nz,nx);
        w_s = zeros(nz,nx);
        u_bd = zeros(nz,nx);
        w_bd = zeros(nz,nx);

        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_s(i,:) = u_st(vec) +  u_st_2(vec);
                w_s(i,:) = w_st(vec) + w_st_2(vec);
                u_bd(i,:) = u_bdt(vec) + u_bdt_2(vec);
                w_bd(i,:) = w_bdt(vec) + w_bdt_2(vec);
        end


        % TE contribution  
        [~,~,u_TEdt,w_TEdt] = DubSorV(muTE(1,i_t),0*muTE(1,i_t),X,Z,[xp(1,i_t) xw(1,i_t)],[zp(1,i_t) zw(1,i_t)],[vtTE(1,i_t,1) vtTE(1,i_t,2)]',[vnTE(1,i_t,1) vnTE(1,i_t,2)]',epSC,SC);

        u_TEdt_2 = 0*u_TEdt;
        w_TEdt_2 = 0*w_TEdt;

        if grd == 1
            [~,~,u_TEdt_2,w_TEdt_2] = DubSorV(muTE(1,i_t),0*muTE(1,i_t),X,Z,[xp(1,i_t) xw(1,i_t)],[-zp(1,i_t) -zw(1,i_t)],[vtTE(1,i_t,1) -vtTE(1,i_t,2)]',[vnTE(1,i_t,1) -vnTE(1,i_t,2)]',epSC,SC);
        end

        u_TEd = zeros(nz,nx);
        w_TEd = zeros(nz,nx);
        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_TEd(i,:) = u_TEdt(vec) + u_TEdt_2(vec);
                w_TEd(i,:) = w_TEdt(vec) + w_TEdt_2(vec);
        end

        % Wake contribution
        [~,~,u_wdt,w_wdt] = DubSorV(muW(1:wakeInd(i_t),i_t)',0*muW(1:wakeInd(i_t),i_t)',X,Z,xw(1:wakeInd(i_t)+1,i_t)',zw(1:wakeInd(i_t)+1,i_t)',[vtw(1:wakeInd(i_t),i_t,1) vtw(1:wakeInd(i_t),i_t,2)]',[vnw(1:wakeInd(i_t),i_t,1) vnw(1:wakeInd(i_t),i_t,2)]',epSC,SC);

        u_wdt_2 = 0*u_wdt;
        w_wdt_2 = 0*w_wdt;

        if grd == 1
            [~,~,u_wdt_2,w_wdt_2] = DubSorV(muW(1:wakeInd(i_t),i_t)',0*muW(1:wakeInd(i_t),i_t)',X,Z,xw(1:wakeInd(i_t)+1,i_t)',-zw(1:wakeInd(i_t)+1,i_t)',[vtw(1:wakeInd(i_t),i_t,1) -vtw(1:wakeInd(i_t),i_t,2)]',[vnw(1:wakeInd(i_t),i_t,1) -vnw(1:wakeInd(i_t),i_t,2)]',epSC,SC);
        end

        u_wd = zeros(nz,nx);
        w_wd = zeros(nz,nx);
        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_wd(i,:) = u_wdt(vec) + u_wdt_2(vec);
                w_wd(i,:) = w_wdt(vec) + w_wdt_2(vec);
        end

        % Lumped vortex contribution
        [~,~,u_lwdt,w_lwdt] = DubSorV(muLump(i_t),0*muLump(i_t),X,Z,[xw(end,i_t) xl(i_t)],[zw(end,i_t) zl(i_t)],[vtlw(1,i_t,1) vtlw(1,i_t,2)]',[vnlw(1,i_t,1) vnlw(1,i_t,2)]',epSC,SC);

        u_lwdt_2 = 0*u_lwdt;
        w_lwdt_2 = 0*w_lwdt;

        if grd == 1
            [~,~,u_lwdt_2,w_lwdt_2] = DubSorV(muLump(i_t),0*muLump(i_t),X,Z,[xw(end,i_t) xl(i_t)],-[zw(end,i_t) zl(i_t)],[vtlw(1,i_t,1) -vtlw(1,i_t,2)]',[vnlw(1,i_t,1) -vnlw(1,i_t,2)]',epSC,SC);
        end

        u_lwd = zeros(nz,nx);
        w_lwd = zeros(nz,nx);
        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_lwd(i,:) = u_lwdt(vec) + u_lwdt_2(vec);
                w_lwd(i,:) = w_lwdt(vec) + w_lwdt_2(vec);
        end

        if LES == 1
            vec = i_t:-1:2;
            u_LEtdt = zeros(1,length(X),length(vec));
            u_LEbdt = zeros(1,length(X),length(vec));
            w_LEtdt = zeros(1,length(X),length(vec));
            w_LEbdt = zeros(1,length(X),length(vec));

            for j = vec
                % LE sheet top contribution
                [~,~,u_LEtdt(1,:,j-1),w_LEtdt(1,:,j-1)] = DubSorV(muLEt(Ncyc*Nstep + 2 - j),0*muLEt(Ncyc*Nstep + 2 - j),X,Z,xpt_LES(:,Ncyc*Nstep + 2 - j)',zpt_LES(:,Ncyc*Nstep + 2 - j),vtLEt(Ncyc*Nstep + 2 - j,:)',vnLEt(Ncyc*Nstep + 2 - j,:)',epSC,SC);

                % LE sheet bot contribution
                [~,~,u_LEbdt(1,:,j-1),w_LEbdt(1,:,j-1)] = DubSorV(muLEb(Ncyc*Nstep + 2 - j),0*muLEb(Ncyc*Nstep + 2 - j),X,Z,xpb_LES(:,Ncyc*Nstep + 2 - j)',zpb_LES(:,Ncyc*Nstep + 2 - j),vtLEb(Ncyc*Nstep + 2 - j,:)',vnLEb(Ncyc*Nstep + 2 - j,:)',epSC,SC);
            end

            u_LEtdt = sum(u_LEtdt,3);
            u_LEbdt = sum(u_LEbdt,3);
            w_LEtdt = sum(w_LEtdt,3);
            w_LEbdt = sum(w_LEbdt,3);

            u_LEtd = zeros(nz,nx);
            w_LEtd = zeros(nz,nx);
            u_LEbd = zeros(nz,nx);
            w_LEbd = zeros(nz,nx);
            for i = 1:nz
                    vec = (i-1)*nx + 1:i*nx;

                    u_LEtd(i,:) = u_LEtdt(vec);
                    w_LEtd(i,:) = w_LEtdt(vec);
                    u_LEbd(i,:) = u_LEbdt(vec);
                    w_LEbd(i,:) = w_LEbdt(vec);
            end
        end

        if LES == 1
            u_d = u_bd + u_TEd + u_wd + u_lwd + u_LEtd + u_LEbd;
            w_d = w_bd + w_TEd + w_wd + w_lwd + w_LEtd + w_LEbd;
        else
            u_d = u_bd + u_TEd + u_wd + u_lwd;
            w_d = w_bd + w_TEd + w_wd + w_lwd;
        end

        u_p = u_d + u_s;
        w_p = w_d + w_s;

        Ut(:,:,figInd) = U + u_p;
        Wt(:,:,figInd) = W + w_p;

        xfi = linspace(-c/4 + x_b(i_t),6*c + x_b(i_t),5*nx)';
        zfi = linspace(-3*c/2,3*c/2,5*nz)';   
       
        [Xfi,Zfi] = meshgrid(xfi,zfi);
        
        u_pi = interp2(Xf,Zf,u_p,Xfi,Zfi);
        w_pi = interp2(Xf,Zf,w_p,Xfi,Zfi);
        
        omega_y = (u_pi(3:end,2:end-1) - u_pi(1:end-2,2:end-1))./(Zfi(3:end,2:end-1) - Zfi(1:end-2,2:end-1)) - (w_pi(2:end-1,3:end) - w_pi(2:end-1,1:end-2))./(Xfi(2:end-1,3:end) - Xfi(2:end-1,1:end-2));
        Xstar = (Xfi(2:end-1,3:end) + Xfi(2:end-1,1:end-2))/2;
        Zstar = (Zfi(3:end,2:end-1) + Zfi(1:end-2,2:end-1))/2;

        if figInd > 1
            close(flowfig)
        end

        

        % Plotting airfoil with LE vortex sheets and the TE vortex
        % sheet.
        flowfig = figure;
        FontSizeAx = 24;
        afFigurePosition = [15 7 25 15];


        set(flowfig, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
        set(flowfig,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

        hold on
        axis equal

        if grd == 1
           plot([min(xp)-10*c; min(xp)+10*c],[0 0],'-k','linewidth',1)
        end

        bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
        redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
        vec = [bluevec; redvec(2:end,:)];
        colormap(vec)

        pcolor(Xstar,Zstar,-omega_y)
%         quiver(Xf,Zf,u_p,w_p,'k')
        plot(xp(:,i_t),zp(:,i_t),'-k','linewidth',2)


        shading interp



        val = 1/10;



%         colorbar
        caxis(5*[-1 1])
        colorbar
        axis([-c/4 + x_b(i_t) 6*c + x_b(i_t) -c c])
%         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
%         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

%         VortMov(i_t-1) = getframe;
        if figInd < 10 
            print('-dpng','-r600',[folder,'/Vorticity/Test_f',num2str(round(f*1000)),'_Wake_00',num2str(figInd)]);
        elseif figInd < 100
            print('-dpng','-r600',[folder,'/Vorticity/Test_f',num2str(round(f*1000)),'_Wake_0',num2str(figInd)]);
        else
            print('-dpng','-r600',[folder,'/Vorticity/Test_f',num2str(round(f*1000)),'_Wake_',num2str(figInd)]);
        end
%         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);

%         print('-dpng','-r300',['FixedV/Vorticity/Valid_St',num2str(St*1000),'_d_c',num2str(d_c*100),'_A_c',num2str(A_c),'_',num2str(figInd),'.png']);
       figInd = figInd + 1;    
    end   
end
save([folder,'/Flowfield',savefilename,'.mat'],'-v7.3','Ut','Wt','Xf','Zf') 

