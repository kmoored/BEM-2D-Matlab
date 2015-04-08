function FlowFieldProcess(f,A_c,d_c,Qinf)

%savefilename = ['_PitchGrd_f',num2str(f),...
%            '_A_c',num2str(A_c),...
%            '_d_c',num2str(d_c)];
%        
%load(['FlowfieldData/Figure9/Processed_',savefilename,'.mat'])  


% Free-stream velocity
alpha = 0;
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

%% Flowfield calculation
Flowfield = 1;
figInd = 1;
epSC = 1e-3;

% Calculating the flowfield around the wing.
nx = 91;
nz = 91;
Ut = zeros(nz,nx,Nstep);
Wt = zeros(nz,nx,Nstep);

if Flowfield == 1
%     for i_t = Ncyc*Nstep+1:Ncyc*Nstep+1
   for i_t = (Ncyc-2)*Nstep+1:Ncyc*Nstep+1
        
        wakeInd = min([(i_t-1) Nlump*Nstep]);
        
        U = Uinf*ones(nz,nx);
        W = Winf*ones(nz,nx);

        % Flow field for ground effect calculations
        xf = linspace(-c/2 + x_b(1,i_t),6*c + x_b(1,i_t),nx)';
        zf = linspace(z_b(1,i_t) - 1.5*c,z_b(1,i_t) + 1.5*c,nz)';

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
        [~,~,u_wdt,w_wdt] = DubSorV(muW(1:wakeInd,i_t)',0*muW(1:wakeInd,i_t)',X,Z,xw(1:wakeInd+1,i_t)',zw(1:wakeInd+1,i_t)',[vtw(1:wakeInd,i_t,1) vtw(1:wakeInd,i_t,2)]',[vnw(1:wakeInd,i_t,1) vnw(1:wakeInd,i_t,2)]',epSC,SC);

        u_wdt_2 = 0*u_wdt;
        w_wdt_2 = 0*w_wdt;

        if grd == 1
            [~,~,u_wdt_2,w_wdt_2] = DubSorV(muW(1:wakeInd,i_t)',0*muW(1:wakeInd,i_t)',X,Z,xw(1:wakeInd+1,i_t)',-zw(1:wakeInd+1,i_t)',[vtw(1:wakeInd,i_t,1) -vtw(1:wakeInd,i_t,2)]',[vnw(1:wakeInd,i_t,1) -vnw(1:wakeInd,i_t,2)]',epSC,SC);
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

        [omega_y,Lambda2,Xstar,Zstar] = Lambda2Crit2D(Xf,Zf,u_p,w_p);

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

        pcolor(Xstar,Zstar,-Lambda2.*omega_y)
%         pcolor(Xstar,Zstar,-omega_y)
        quiver(Xf,Zf,u_p,w_p,'k')
        plot(xp(:,i_t),zp(:,i_t),'-k','linewidth',2)

%         if grd == 1
%             plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
%         end

        shading interp

        if LES == 1
            plot(xplus_t,zplus_t,'xk')
            plot(xminus_t,zminus_t,'xk')
            plot(xplus_b,zplus_b,'xk')
            plot(xminus_b,zminus_b,'xk')
            quiver(xp(Indsep_t,i_t),zp(Indsep_t,i_t),Vshear_t(1,i_t+1)*delT,Vshear_t(2,i_t+1)*delT,'r')
            quiver(xp(Indsep_b,i_t),zp(Indsep_b,i_t),Vshear_b(1,i_t+1)*delT,Vshear_b(2,i_t+1)*delT,'r')
        end
%         plot(xTE(:,i_t),zTE(:,i_t),'.-b','linewidth',2)

        val = 1/10;

        if LES == 1 && t > 0
            plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
            plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)

            for j = i_t:-1:2
                if muLEt(Ncyc*Nstep + 2 - j) == 0
                else
                    plot(xpt_LES(:,Ncyc*Nstep + 2 - j),zpt_LES(:,Ncyc*Nstep + 2 - j),'.-g','linewidth',1)
    %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vtLEt(Ncyc*Nstep + 1 - j,1),val*vtLEt(Ncyc*Nstep + 1 - j,2),'b')
    %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vnLEt(Ncyc*Nstep + 1 - j,1),val*vnLEt(Ncyc*Nstep + 1 - j,2),'k')

                end
            end

            for j = i_t:-1:2
                if muLEb(Ncyc*Nstep + 2 - j) == 0   
                else
                    plot(xpb_LES(:,Ncyc*Nstep + 2 - j),zpb_LES(:,Ncyc*Nstep + 2 - j),'.-r','linewidth',1)
    %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vtLEb(Ncyc*Nstep + 1 - j,1),val*vtLEb(Ncyc*Nstep + 1 - j,2),'b')
    %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vnLEb(Ncyc*Nstep + 1 - j,1),val*vnLEb(Ncyc*Nstep + 1 - j,2),'k')

                end
            end

        end

        if i_t > 1
%             plot(xp(Stagpt,i_t),zp(Stagpt,i_t),'xk','linewidth',2)
%             plot(xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),'.-b')
        end



%         colorbar
        caxis(6*[-1 1])
        colorbar
        axis([-c/2 + x_b(1,i_t) 6*c + x_b(1,i_t) z_b(1,i_t) - 1.5*c z_b(1,i_t) + 1.5*c])
%         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
%         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

%         VortMov(i_t-1) = getframe;
%         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);
% 
%         if figInd < 10 
%             print('-depsc','-r600',['FixedV/A_c_',num2str(A_c),'/St_',num2str(St),'/Vort/Vort_00',num2str(figInd)]);
%         elseif figInd < 100
%             print('-depsc','-r600',['FixedV/A_c_',num2str(A_c),'/St_',num2str(St),'/Vort/Vort_0',num2str(figInd)]);
%         else
%             print('-depsc','-r600',['FixedV/A_c_',num2str(A_c),'/St_',num2str(St),'/Vort/Vort_',num2str(figInd)]);
%         end
       figInd = figInd + 1;    
    end   
end
save(['FlowfieldData/Figure9/Flowfield',savefilename,'.mat'],'-v7.3','Ut','Wt')    

% %% Plotting time-averaged flowfield
% figure;
% FontSizeAx = 24;
% afFigurePosition = [15 7 25 15];
% 
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% axis equal
% 
% % if grd == 1
% %    plot([min(xp)-10*c; min(xp)+10*c],[0 0],'-k','linewidth',1)
% % end
% 
% bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
% redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
% vec = [bluevec; redvec(2:end,:)];
% colormap(vec)
% 
% quiver(Xf,Zf,mean(Ut-Uinf,3),mean(Wt-Winf,3),'k')
% % quiver(Xf,Zf,Uinf*ones(nz,nx),mean(Wt(:,:,1:151)-Winf,3),'k')
% plot(xp_0 ,zp_0 + d_c*c,'-k','linewidth',2)
% hold off
% 
% if grd == 1
%     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% end
% 
% shading interp
% 
% 
% 
% %         colorbar
% caxis(6*[-1 1])
% colorbar
% axis([-c/4 3*c  0 2*d_c*c])
% %         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
% %         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% %         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% 
% %         VortMov(i_t-1) = getframe;
% %         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);
% 
%             print('-dpng','-r300',['D:/Data/grd_2D/v4/Up/TAFF',num2str(St*1000),'_d_c',num2str(d_c*100),'_A_c',num2str(A_c),'_',num2str(figInd),'.png']);
 
