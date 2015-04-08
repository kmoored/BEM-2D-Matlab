clear
close all
clc

% Parameters for identifying data file.
St = 0.25;
A_c = 0.05;
Qinf = 4700;

savefilename = ['_Pitch_St',num2str(St),...
        '_A_c',num2str(A_c),...
        '_Re',num2str(Re)];

%% Loading data 
load(['Theodorson/Parameters',savefilename,'.mat']);
Data = load(['Theodorson/Data',savefilename,'.txt']);
PanelProp = load(['Theodorson/PanelProp',savefilename,'.txt']);
WakeProp = load(['Theodorson/WakeProp',savefilename,'.txt']);

%% Defining variables

% Free-stream velocity
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

wakeInd = Data(:,1)';
Q0(1,:) = Data(:,2)';
Q0(2,:) = Data(:,3)';
x_b = Data(:,4)';
z_b = Data(:,5)';
a_b = Data(:,6)';
Fx = Data(:,7)';
Fz = Data(:,8)';
Pow = Data(:,9)';
L = Data(:,10)';
T = Data(:,11)';
D_visc = Data(:,12)';
Gamma = Data(:,13)';
Cf = Data(:,14)';
Cl = Data(:,15)';
Ct = Data(:,16)';
Cpow = Data(:,17)';
Fx_s = Data(:,18)';
Fx_us = Data(:,19)';
Fz_s = Data(:,20)';
Fz_us = Data(:,21)';
Pow_s = Data(:,22)';
Pow_us = Data(:,23)';
L_s = Data(:,24)';
L_us = Data(:,25)';
T_s = Data(:,26)';
T_us = Data(:,27)';
Cl_s = Data(:,28)';
Cl_us = Data(:,29)';
Ct_s = Data(:,30)';
Ct_us = Data(:,31)';
Cpow_s = Data(:,32)';
Cpow_us = Data(:,33)';

clear Data

xp_0 = reshape(PanelProp(:,1),[Npanels (Ncyc*Nstep + 1)]);
xp_0 = [xp_0;xp_0(1,:)];
zp_0 = reshape(PanelProp(:,2),[Npanels (Ncyc*Nstep + 1)]);
zp_0 = [zp_0;zp_0(1,:)];
xp = reshape(PanelProp(:,3),[Npanels (Ncyc*Nstep + 1)]);
xp = [xp;xp(1,:)];
zp = reshape(PanelProp(:,4),[Npanels (Ncyc*Nstep + 1)]);
zp = [zp;zp(1,:)];
Vc(:,:,1)  = reshape(PanelProp(:,5),[Npanels (Ncyc*Nstep + 1)]);
Vc(:,:,2) = reshape(PanelProp(:,6),[Npanels (Ncyc*Nstep + 1)]);
xc = reshape(PanelProp(:,7),[Npanels (Ncyc*Nstep + 1)]);
zc = reshape(PanelProp(:,8),[Npanels (Ncyc*Nstep + 1)]);
vt(:,:,1) = reshape(PanelProp(:,9),[Npanels (Ncyc*Nstep + 1)]);
vt(:,:,2) = reshape(PanelProp(:,10),[Npanels (Ncyc*Nstep + 1)]);
vn(:,:,1) = reshape(PanelProp(:,11),[Npanels (Ncyc*Nstep + 1)]);
vn(:,:,2) = reshape(PanelProp(:,12),[Npanels (Ncyc*Nstep + 1)]);
dL = reshape(PanelProp(:,13),[Npanels (Ncyc*Nstep + 1)]);
Xc = reshape(PanelProp(:,14),[Npanels (Ncyc*Nstep + 1)]);
Zc = reshape(PanelProp(:,15),[Npanels (Ncyc*Nstep + 1)]);
sigma = reshape(PanelProp(:,16),[Npanels (Ncyc*Nstep + 1)]);
mu = reshape(PanelProp(:,17),[Npanels (Ncyc*Nstep + 1)]);
Qp = reshape(PanelProp(:,18),[Npanels (Ncyc*Nstep + 1)]);
Qt = reshape(PanelProp(:,19),[Npanels (Ncyc*Nstep + 1)]);
Cp_s = reshape(PanelProp(:,20),[Npanels (Ncyc*Nstep + 1)]);
Cp_us = reshape(PanelProp(:,21),[Npanels (Ncyc*Nstep + 1)]);
Cp = reshape(PanelProp(:,22),[Npanels (Ncyc*Nstep + 1)]);
dFshear = reshape(PanelProp(:,23),[Npanels (Ncyc*Nstep + 1)]);
      
clear PanelProp

xwtemp = reshape(WakeProp(:,1),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
xw = xwtemp(1:end-1,:);
xl = xwtemp(end,:);
zwtemp = reshape(WakeProp(:,2),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
zw = zwtemp(1:end-1,:);
zl = zwtemp(end,:);
vtwtemp_1 = reshape(WakeProp(:,3),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vtwtemp_2 = reshape(WakeProp(:,4),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vtTE(1,:,1) = vtwtemp_1(1,:);
vtTE(1,:,2) = vtwtemp_2(1,:);
vtw(:,:,1) = vtwtemp_1(2:end-1,:);
vtw(:,:,2) = vtwtemp_2(2:end-1,:);
vtlw(1,:,1) = vtwtemp_1(end,:);
vtlw(1,:,2) = vtwtemp_2(end,:);
vnwtemp_1 = reshape(WakeProp(:,5),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vnwtemp_2 = reshape(WakeProp(:,6),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vnTE(1,:,1) = vnwtemp_1(1,:);
vnTE(1,:,2) = vnwtemp_2(1,:);
vnw(:,:,1) = vnwtemp_1(2:end-1,:);
vnw(:,:,2) = vnwtemp_2(2:end-1,:);
vnlw(1,:,1) = vnwtemp_1(end,:);
vnlw(1,:,2) = vnwtemp_2(end,:);
muWtemp = reshape(WakeProp(:,7),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
muTE(1,:) = muWtemp(1,:);
muW = muWtemp(2:end-1,:);
muLump(1,:) = muWtemp(end,:);
GammaWtemp = reshape(WakeProp(:,8),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
GammaW = GammaWtemp(1:end-1,:);
GammaLump = GammaWtemp(end,:);

clear WakeProp xwtemp zwtemp vtwtemp_1 vtwtemp_2 vnwtemp_1 vnwtemp_2 muWtemp GammaWtemp


%% Flowfield calculation
Flowfield = 1;
figInd = 1;
epSC = 2.5e-3;

% Calculating the flowfield around the wing.
nx = 51;
nz = 51;
Ut = zeros(nz,nx,Nstep);
Wt = zeros(nz,nx,Nstep);

if Flowfield == 1
    for i_t = (Ncyc-1)*Nstep+1:Ncyc*Nstep+1
        U = Uinf*ones(nz,nx);
        W = Winf*ones(nz,nx);

        % Flow field for ground effect calculations
        xf = linspace(-c/4 + x_b(i_t),4*c + x_b(i_t),nx)';
        zf = linspace(0,3*c/2,nz)';

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

        omega_y = (u_p(3:end,2:end-1) - u_p(1:end-2,2:end-1))./(Zf(3:end,2:end-1) - Zf(1:end-2,2:end-1)) - (w_p(2:end-1,3:end) - w_p(2:end-1,1:end-2))./(Xf(2:end-1,3:end) - Xf(2:end-1,1:end-2));
        Xstar = (Xf(2:end-1,3:end) + Xf(2:end-1,1:end-2))/2;
        Zstar = (Zf(3:end,2:end-1) + Zf(1:end-2,2:end-1))/2;

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
        axis([-c/4 + x_b(i_t) 4*c + x_b(i_t) 0 c])
%         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
%         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

%         VortMov(i_t-1) = getframe;
%         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);

        print('-dpng','-r300',['FixedV/Vorticity/Valid_St',num2str(St*1000),'_d_c',num2str(d_c*100),'_A_c',num2str(A_c),'_',num2str(figInd),'.png']);
       figInd = figInd + 1;    
    end   
end
% save(['GrdData/Flowfield',savefilename,'.mat'],'-v7.3','Ut','Wt')    

%% Plotting time-averaged flowfield
figure;
FontSizeAx = 24;
afFigurePosition = [15 7 25 15];

set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
axis equal

if grd == 1
   plot([min(xp)-10*c; min(xp)+10*c],[0 0],'-k','linewidth',1)
end

bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
vec = [bluevec; redvec(2:end,:)];
colormap(vec)

quiver(Xf,Zf,mean(Ut-Uinf,3),mean(Wt-Winf,3),'k')
% quiver(Xf,Zf,Uinf*ones(nz,nx),mean(Wt(:,:,1:151)-Winf,3),'k')
plot(xp_0 + x_b(end) ,zp_0 + d_c*c,'-k','linewidth',2)
hold off

%         if grd == 1
%             plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
%         end

shading interp



%         colorbar
caxis(6*[-1 1])
colorbar
axis([-c/4 + x_b(i_t) 3*c + x_b(i_t) 0 c])
%         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
%         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

%         VortMov(i_t-1) = getframe;
%         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);

%             print('-dpng','-r300',['GrdData/Vorticity/GrdEffect_St',num2str(St*1000),'_d_c',num2str(d_c*100),'_A_c',num2str(A_c),'_',num2str(figInd),'.png']);
        
%% Plotting the forces acting on the body
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [0 2 30 20];
axespos = [0.175 0.15 0.775 0.8];
ylabelpos = [-0.125 0.5];
xlabelpos = [0.5 -0.1];

t = ([1:Ncyc*Nstep+1] - 1)*delT;  
delt = (t((Ncyc-2)*Nstep+1:end) - t((Ncyc-2)*Nstep+1))*f;
span = 0.05;

% Lift
figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(delt,Cl_s((Ncyc-2)*Nstep+1:end),'-b','linewidth',3)
    plot(delt,smooth(Cl_us((Ncyc-2)*Nstep+1:end),span,'loess'),'-r','linewidth',3)
    plot(delt,smooth(Cl((Ncyc-2)*Nstep+1:end),span,'loess'),'-k','linewidth',3)

    plot([min(delt) max(delt)],mean(Cl((Ncyc-2)*Nstep+1:end))*[1 1],':k','linewidth',2)
    plot([min(delt) max(delt)],[0 0],'-k','linewidth',2)
hold off

axis([0 2 1.2*min(Cl) 1.2*max(Cl)])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Steady','Unsteady','Total','Location','NorthWest')
xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_l$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

% Thrust
figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(delt,Ct_s((Ncyc-2)*Nstep+1:end),'-b','linewidth',3)
    plot(delt,smooth(Ct_us((Ncyc-2)*Nstep+1:end),span,'loess'),'-r','linewidth',3)
    plot(delt,smooth(Ct((Ncyc-2)*Nstep+1:end),span,'loess'),'-k','linewidth',3)
    
    plot([min(delt) max(delt)],mean(Ct((Ncyc-2)*Nstep+1:end))*[1 1],':k','linewidth',2)
    plot([min(delt) max(delt)],[0 0],'-k','linewidth',2)
hold off

axis([0 2 1.2*min(Ct) 1.2*max(Ct)])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Steady','Unsteady','Total','Location','SouthEast')
xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

% Circulation
figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(delt,(-2/c/Uinf)*Gamma((Ncyc-2)*Nstep+1:end),'-k','linewidth',3)
    
    plot([min(delt) max(delt)],mean((-2/c/Uinf)*Gamma((Ncyc-2)*Nstep+1:end))*[1 1],':k','linewidth',2)
    plot([min(delt) max(delt)],mean(Cl_s((Ncyc-2)*Nstep+1:end))*[1 1],'--b','linewidth',2)
    plot([min(delt) max(delt)],mean(Cl_us((Ncyc-2)*Nstep+1:end))*[1 1],'--g','linewidth',2)
    plot([min(delt) max(delt)],[0 0],'-k','linewidth',2)
hold off

axis([0 2 1.2*min((-2/c/Uinf)*Gamma) 1.2*max((-2/c/Uinf)*Gamma)])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Circulatory force','Time-avg circulatory force','Time-avg steady force','Time-avg unsteady force','Location','SouthEast')
xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$-\frac{2\Gamma}{cU_{\infty}}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');


% Pressure coefficient
figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

xc_0 = (xp_0(1:end-1,1)/2 + xp_0(2:end,1)/2);
hold on
    plot(xc_0(1:Npanels/2)/c,mean(Cp(1:Npanels/2,(Ncyc-2)*Nstep+1:end),2),'ob','linewidth',3,'markersize',6)
    plot(xc_0(Npanels/2+1:end)/c,mean(Cp(Npanels/2+1:end,(Ncyc-2)*Nstep+1:end),2),'sk','linewidth',3,'markersize',6)
hold off

axis([0 1 1.2*min(mean(Cp(:,(Ncyc-2)*Nstep+1:end),2)) 1.2*max(mean(Cp(:,(Ncyc-2)*Nstep+1:end),2))])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Bottom','Top','Location','South')
xlabel('$$x/c$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\bar{C_p}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');



%% Airfoil plotting
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [0 2 30 20];
% axespos = [0.15 0.2 0.84 0.8];
% ylabelpos = [-0.12 0.5];
% xlabelpos = [0.5 -0.25];


% % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
% fighand = figure;
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% axis equal
% 
% plot(Xc(:,end),Zc(:,end),'xr','linewidth',2)
% plot(Xc(1,end),Zc(1,end),'xb','linewidth',2)
% plot(Xc(end,end),Zc(end,end),'xg','linewidth',2)
% plot(xp(:,end),zp(:,end),'.-k','linewidth',2)
% plot(xTE(:,end),zTE(:,end),'.-b','linewidth',2)
% % quiver(xc(:,end),zc(:,end),vn(:,1,end),vn(:,2,end),'b','linewidth',2)
% val = 1/10;
% 
% if grd == 1
% %         plot(xp_2(:,i_t),zp_2(:,i_t),'-k','linewidth',2)
% %         plot(xTE_2(:,i_t),zTE_2(:,i_t),'.-b','linewidth',2)
%     plot([min(xp(:,end))-10*c; min(xp(:,end))+10*c],[0 0],'-k','linewidth',4)
% %         if t > 0
% %             plot(xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),'.-b','linewidth',2,'markersize',14)
% %         end
% end
% 
% if t > 0
%     if LES == 1
%         plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
%         plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)
% 
%         for j = i_t:-1:2
%             if muLEt(Ncyc*Nstep + 2 - j) == 0
%             else
%                 plot(xpt_LES(:,Ncyc*Nstep + 2 - j),zpt_LES(:,Ncyc*Nstep + 2 - j),'.-g','linewidth',2,'markersize',14)
% %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vtLEt(Ncyc*Nstep + 1 - j,1),val*vtLEt(Ncyc*Nstep + 1 - j,2),'b')
% %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vnLEt(Ncyc*Nstep + 1 - j,1),val*vnLEt(Ncyc*Nstep + 1 - j,2),'k')
% 
%             end
%         end
% 
%         for j = i_t:-1:2
%             if muLEb(Ncyc*Nstep + 2 - j) == 0  
%             else
%                 plot(xpb_LES(:,Ncyc*Nstep + 2 - j),zpb_LES(:,Ncyc*Nstep + 2 - j),'.-r','linewidth',2,'markersize',14)
% %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vtLEb(Ncyc*Nstep + 1 - j,1),val*vtLEb(Ncyc*Nstep + 1 - j,2),'b')
% %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vnLEb(Ncyc*Nstep + 1 - j,1),val*vnLEb(Ncyc*Nstep + 1 - j,2),'k')
% 
%             end
%         end
%     end
% %         plot(xp(Stagpt,i_t),zp(Stagpt,i_t),'xk','linewidth',2)
% 
%     bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
%     redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
%     vec = [bluevec; redvec(2:end,:)];
%     vec = vec(end:-1:1,:);
% 
%     WakeCirc = GammaW'/max(GammaW)/(1/2);
%     WakeCirc(WakeCirc > 1) = 1;  
%     WakeCirc(WakeCirc < -1) = -1;
%     WakeCirc = WakeCirc + 1;
%     WakeCirc = round(WakeCirc*100) + 1;
%     WakeCirc(WakeCirc > 199) = 199;
% 
%     for i_w = 1:Ncyc*Nstep+1
%         plot(xw(i_w),zw(i_w),'.','color',vec(WakeCirc(i_w),:),'linewidth',2,'markersize',14)
%     end
% end
% 
% axis equal
% %     axis([-c/2 + x_b(i_t) 2*c + x_b(i_t) -1/2*c + z_b(i_t) 1/2*c + z_b(i_t)])
% axis([-c/2 + x_b(end-1) 7*c + x_b(end-1) z_b(end-1) - 1*c z_b(end-1) + 1*c])
% 
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$x$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$z$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');



%% Analytical solution.

% [~,~,~,~,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,600,tmax/2);
% % 
% % xp_a = [xpb;xpt(2:end)];
% % % zp = [zpb;zpt(2:end)];
% % zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];
% 
% Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
% Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
% D0 = a_vdv*(1 - k + k*epsilon);
% D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
% D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);
% 
% u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
% w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));
% 
% Cp_anltc_t = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;
% 
% alpha_b = -alpha;
% [xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,600,tmax/2);
% 
% xp_a = [xpb;xpt(2:end)];
% % zp = [zpb;zpt(2:end)];
% zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];
% 
% Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
% Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
% D0 = a_vdv*(1 - k + k*epsilon);
% D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
% D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);
% 
% u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
% w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));
% 
% Cp_anltc_b = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;


% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp(:,end),zp(:,end),'.-k')
% plot(xTE(:,end),zTE(:,end),'.-g')
% plot(xw,zw,'.-b')
% if grd == 1
%     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% end
% axis equal
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% axis([-c/2 + x_b(i_t) 4*c + x_b(i_t) -c/10 c])


% figure
% hold on
% plot(((1:i_t) - 1)*delT,x_b(1:end-1),'-b','linewidth',2)
% xlabel('$$t$$, sec','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% ylabel('Position','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
      

% % Plotting pressure coefficient both numerical and analytical solutions.
% figure
% hold on
% 
% 
% plot((xc(1:Npanels/2,end) - x_b(end-1))/c,-Cp(1:Npanels/2,end-round(Nstep/4)),'sk')
% plot((xc(Npanels/2 + 1:end,end) - x_b(end-1))/c,-Cp(Npanels/2 + 1:end,end-round(Nstep/4)),'^k')
% 
% % plot(xp_a(round(599/2):end)/c,-Cp_anltc_t,'-k')
% % plot(xp_a(round(599/2):end)/c,-Cp_anltc_b,'-k')
% 
% % plot(xp(2:round((Npanels + 1)/2)),-Cp(1:round((Npanels+1)/2)-1),'sk')
% % plot(xp(round((Npanels + 1)/2):end-1),-Cp(round((Npanels+1)/2)-1:end),'^k')
% xlabel('$$\frac{x}{c}$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% ylabel('$$-C_p$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% % axis([0 c -1 2])
% % axis([-1.6 -1.3 -1 2])
% legend('Bottom','Top','Location','NorthEast')
% 
% 


%% Plotting coefficient of pressure.
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [1 1 20 15];
% ylabelpos = [-0.1 0.5];
% xlabelpos = [0.5 -0.15];
% axespos = [0.15 0.2 0.8 0.75];
% Panelt = 189;
% Panelb = 21;
% 
% t = ([1:Ncyc*Nstep+1] - 1)*delT;  
% delt = (t((Ncyc-1)*Nstep+1:end) - t((Ncyc-1)*Nstep+1))*f;
% 
% Press = (Cp(Panelt,(Ncyc-1)*Nstep+1:end) + 1)*(1/2*rho*Qinf^2);
% Press_avg = mean(Press);
% Press_var = Press - Press_avg;
% 
% Press_max_t = (max(Cp(Npanels/2 + 1:end,(Ncyc-1)*Nstep+1:end),[],2))*(1/2*rho*Qinf^2);
% Press_max = Press_max_t - mean(Press_max_t);
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% %     plot(delt,Cp_s(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--b','linewidth',3)
% %     plot(delt,Cp_us(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--r','linewidth',3)
%     plot(delt,Press_var,'-k','linewidth',3)
% hold off
% 
% axis([0 1 1.2*min(Press_var) 1.2*max(Press_var)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$P$$ (Pa)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% 
% 
% 
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% %     plot(delt,Cp_s(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--b','linewidth',3)
% %     plot(delt,Cp_us(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--r','linewidth',3)
%     plot(xc(Npanels/2+1:end,1)/c,Press_max,'-k','linewidth',3)
% hold off
% 
% axis([0 1 1.2*min(Press_max) 1.2*max(Press_max)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$x/c$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$P$$ (Pa)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

% 
% 
% 
% 
% % Plotting of the forces acting on the body.
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(delt,L((Ncyc-1)*Nstep+1:end),'-r','linewidth',3)
%     plot(delt,T((Ncyc-1)*Nstep+1:end),'-b','linewidth',3)
% hold off
% 
% % axis([0 1 1.2*min(Press_var) 1.2*max(Press_var)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$F$$ (N)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% 
% % Calculating the moment from each panel about the leading edge.
% % Leading edge origin:
% x0 = kron(xp(Npanels/2 + 1,:),ones(Npanels,1));
% z0 = kron(zp(Npanels/2 + 1,:),ones(Npanels,1));
% r0(:,1,:) = x0;
% r0(:,2,:) = 0*z0;
% r0(:,3,:) = z0;
% r(:,1,:) = xc;
% r(:,2,:) = 0*zc;
% r(:,3,:) = zc;
% delF3(:,1,:) = delF(:,1,:);
% delF3(:,2,:) = 0*delF(:,1,:);
% delF3(:,3,:) = delF(:,2,:);
% 
% dM = cross(r-r0,delF3);
% dM = dM(:,2,:);
% Torque(1,:) = -sum(dM,1);
% theta_dot = -2*pi*f*alpha_max*cos(2*pi*f*t);
% Pow_T = abs(Torque.*theta_dot);
% Pow_T_avg = mean(Pow_T((Ncyc-1)*Nstep+1:end))
% 
% C_pow_T = Pow_T_avg/(1/2*rho*Qinf^3*c*b)
% 
% % Plotting of the forces acting on the body.
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [1 1 20 15];
% ylabelpos = [-0.125 0.5];
% xlabelpos = [0.5 -0.15];
% axespos = [0.175 0.2 0.8 0.75];
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(delt,Pow((Ncyc-1)*Nstep+1:end),'-r','linewidth',3)
%      plot([0 1],[Pow_avg Pow_avg],'-r','linewidth',1)
%     plot(delt,Pow_T((Ncyc-1)*Nstep+1:end),'-b','linewidth',3)
%     plot([0 1],[Pow_T_avg Pow_T_avg],'-b','linewidth',1)
% hold off
% 
% % axis([0 1 1.2*min(Press_var) 1.2*max(Press_var)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$Pow$$ (W)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

% %% Intantaneous Power
% 
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [0 2 30 20];
% axespos = [0.12 0.15 0.84 0.8];
% ylabelpos = [-0.075 0.5];
% xlabelpos = [0.5 -0.1];
% 
% t = ([1:Ncyc*Nstep+1] - 1)*delT;  
% delt = (t((Ncyc-2)*Nstep+1:end) - t((Ncyc-2)*Nstep+1))*f;
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% a_TE = [(zp(end,(Ncyc-2)*Nstep+1:end) - zp(Npanels/2+1,(Ncyc-2)*Nstep+1:end))/A_TE]';
% 
% Cpow_last = Cpow((Ncyc-2)*Nstep+1:end);
% hold on
%     plot(1/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(2/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(3/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(4/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',2)
%     plot(5/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(6/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(7/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
% 
%     
%     plot([0 2],[0 0],'-k','linewidth',1)
%     plot(delt,3*a_TE,'--k','linewidth',3)
%     plot(delt,Cpow_last,'-b','linewidth',3)
% hold off
% 
% axis([0 2 1.2*min(Cpow_last) 1.2*max(Cpow_last)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

%% Intantaneous Lift
% 
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [0 2 30 20];
% axespos = [0.12 0.15 0.84 0.8];
% ylabelpos = [-0.075 0.5];
% xlabelpos = [0.5 -0.1];
% 
% t = ((1:Ncyc*Nstep+1) - 1)*delT;  
% delt = t*f;
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % a_TE = [(zp(end,(Ncyc-2)*Nstep+1:end) - zp(Npanels/2+1,(Ncyc-2)*Nstep+1:end))/A_TE]';
% 
% hold on   
%     plot([0 Ncyc],[0 0],'-k','linewidth',1)
%     plot(delt,Ct_s,'-b','linewidth',3)
%     plot(delt,Ct_us,'-r','linewidth',3)
%     plot(delt,Ct,'-k','linewidth',3)
% hold off
% 
% axis([0 Ncyc min(Ct_us) max(Ct_us)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_l$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

%% Calculating the flowfield around the wing.
% nx = 50;
% nz = 50;
% 
% U = Uinf*ones(nz,nx);
% W = Winf*ones(nz,nx);
% 
% %         xf = linspace(-c/2 + x_b(i_t),3*c + x_b(i_t),nx)';
% %         zf = linspace(-3*c/4 + z_b(i_t),3*c/4 + z_b(i_t),nz)';
% 
% % Flow field for ground effect calculations
% xf = linspace(-c/2 + x_b(i_t),3*c + x_b(i_t),nx)';
% zf = linspace(0,c,nz)';
% 
% [Xf,Zf] = meshgrid(xf,zf);
% 
% X = zeros(1,nx*nz);
% Z = zeros(1,nx*nz);
% 
% for i = 1:nz
%     vec = (i-1)*nx + 1:i*nx;
% 
%     X(1,vec) = Xf(i,:);
%     Z(1,vec) = Zf(i,:);
% end
% 
% % Calculating flow field
% 
% % Body contribution
% [u_st,w_st,u_bdt,w_bdt] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',epSC,SC);
% 
% u_st_2 = 0*u_st;
% w_st_2 = 0*w_st;
% u_bdt_2 = 0*u_bdt;
% w_bdt_2 = 0*w_bdt;
% 
% if grd == 1
%     [u_st_2,w_st_2,u_bdt_2,w_bdt_2] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp_2(:,i_t)',zp_2(:,i_t)',vt_2(:,:,i_t)',vn_2(:,:,i_t)',epSC,SC);
% end
% 
% u_s = zeros(nz,nx);
% w_s = zeros(nz,nx);
% u_bd = zeros(nz,nx);
% w_bd = zeros(nz,nx);
% 
% for i = 1:nz
%         vec = (i-1)*nx + 1:i*nx;
% 
%         u_s(i,:) = u_st(vec) +  u_st_2(vec);
%         w_s(i,:) = w_st(vec) + w_st_2(vec);
%         u_bd(i,:) = u_bdt(vec) + u_bdt_2(vec);
%         w_bd(i,:) = w_bdt(vec) + w_bdt_2(vec);
% end
% 
% 
% % TE contribution  
% [~,~,u_TEdt,w_TEdt] = DubSorV(muTE(i_t),0*muTE(i_t),X,Z,xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',epSC,SC);
% 
% u_TEdt_2 = 0*u_TEdt;
% w_TEdt_2 = 0*w_TEdt;
% 
% if grd == 1
%     [~,~,u_TEdt_2,w_TEdt_2] = DubSorV(muTE(i_t),0*muTE(i_t),X,Z,xTE_2(:,i_t)',zTE_2(:,i_t)',vtTE_2(1,:,i_t)',vnTE_2(1,:,i_t)',epSC,SC);
% end
% 
% u_TEd = zeros(nz,nx);
% w_TEd = zeros(nz,nx);
% for i = 1:nz
%         vec = (i-1)*nx + 1:i*nx;
% 
%         u_TEd(i,:) = u_TEdt(vec) + u_TEdt_2(vec);
%         w_TEd(i,:) = w_TEdt(vec) + w_TEdt_2(vec);
% end
% 
% % Wake contribution
% [~,~,u_wdt,w_wdt] = DubSorV(muW(Ncyc*Nstep + 2 - i_t:end),0*muW(Ncyc*Nstep + 2 - i_t:end),X,Z,xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',epSC,SC);
% 
% u_wdt_2 = 0*u_wdt;
% w_wdt_2 = 0*w_wdt;
% 
% if grd == 1
%     [~,~,u_wdt_2,w_wdt_2] = DubSorV(muW(Ncyc*Nstep + 2 - i_t:end),0*muW(Ncyc*Nstep + 2 - i_t:end),X,Z,xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),vtw_2(Ncyc*Nstep + 2 - i_t:end,:)',vnw_2(Ncyc*Nstep + 2 - i_t:end,:)',epSC,SC);
% end
% 
% u_wd = zeros(nz,nx);
% w_wd = zeros(nz,nx);
% for i = 1:nz
%         vec = (i-1)*nx + 1:i*nx;
% 
%         u_wd(i,:) = u_wdt(vec) + u_wdt_2(vec);
%         w_wd(i,:) = w_wdt(vec) + w_wdt_2(vec);
% end
% 
% if LES == 1
%     vec = i_t:-1:2;
%     u_LEtdt = zeros(1,length(X),length(vec));
%     u_LEbdt = zeros(1,length(X),length(vec));
%     w_LEtdt = zeros(1,length(X),length(vec));
%     w_LEbdt = zeros(1,length(X),length(vec));
% 
%     for j = vec
%         % LE sheet top contribution
%         [~,~,u_LEtdt(1,:,j-1),w_LEtdt(1,:,j-1)] = DubSorV(muLEt(Ncyc*Nstep + 2 - j),0*muLEt(Ncyc*Nstep + 2 - j),X,Z,xpt_LES(:,Ncyc*Nstep + 2 - j)',zpt_LES(:,Ncyc*Nstep + 2 - j),vtLEt(Ncyc*Nstep + 2 - j,:)',vnLEt(Ncyc*Nstep + 2 - j,:)',epSC,SC);
% 
%         % LE sheet bot contribution
%         [~,~,u_LEbdt(1,:,j-1),w_LEbdt(1,:,j-1)] = DubSorV(muLEb(Ncyc*Nstep + 2 - j),0*muLEb(Ncyc*Nstep + 2 - j),X,Z,xpb_LES(:,Ncyc*Nstep + 2 - j)',zpb_LES(:,Ncyc*Nstep + 2 - j),vtLEb(Ncyc*Nstep + 2 - j,:)',vnLEb(Ncyc*Nstep + 2 - j,:)',epSC,SC);
%     end
% 
%     u_LEtdt = sum(u_LEtdt,3);
%     u_LEbdt = sum(u_LEbdt,3);
%     w_LEtdt = sum(w_LEtdt,3);
%     w_LEbdt = sum(w_LEbdt,3);
% 
%     u_LEtd = zeros(nz,nx);
%     w_LEtd = zeros(nz,nx);
%     u_LEbd = zeros(nz,nx);
%     w_LEbd = zeros(nz,nx);
%     for i = 1:nz
%             vec = (i-1)*nx + 1:i*nx;
% 
%             u_LEtd(i,:) = u_LEtdt(vec);
%             w_LEtd(i,:) = w_LEtdt(vec);
%             u_LEbd(i,:) = u_LEbdt(vec);
%             w_LEbd(i,:) = w_LEbdt(vec);
%     end
% end
% 
% if LES == 1
%     u_d = u_bd + u_TEd + u_wd + u_LEtd + u_LEbd;
%     w_d = w_bd + w_TEd + w_wd + w_LEtd + w_LEbd;
% else
%     u_d = u_bd + u_TEd + u_wd;
%     w_d = w_bd + w_TEd + w_wd;
% end
% 
% 
% u_p = u_d + u_s;
% w_p = w_d + w_s;
% 
% U = U + u_p;
% W = W + w_p;
% 
% omega_y = (u_p(3:end,2:end-1) - u_p(1:end-2,2:end-1))./(Zf(3:end,2:end-1) - Zf(1:end-2,2:end-1)) - (w_p(2:end-1,3:end) - w_p(2:end-1,1:end-2))./(Xf(2:end-1,3:end) - Xf(2:end-1,1:end-2));
% Xstar = (Xf(2:end-1,3:end) + Xf(2:end-1,1:end-2))/2;
% Zstar = (Zf(3:end,2:end-1) + Zf(1:end-2,2:end-1))/2;
% 
% 
% 
% % % Plotting discretized airfoil with normals and lift vectors.
% % figure
% % set(gcf,'DefaultAxesfontsize',20,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% % 
% % hold on
% % plot(xp(:,end),zp(:,end),'.-k','linewidth',2)
% % plot(xTE(:,end),zTE(:,end),'.-g')
% % plot(xw,zw,'.-g')
% % plot(Xc,Zc,'xb')
% % if grd == 1
% %     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% % end
% % streamline(Xf,Zf,U,W,Xf(:,1),Zf(:,1))
% % % quiver(Xf,Zf,U,W)
% % axis equal
% % axis([-c/2 + x_b(end-1) 2*c + x_b(end-1) -c/10 c])
% % % axis([-c/2 + x_b(end-1) 2*c + x_b(end-1) -c/2 + z_b(end-1) c/2 + z_b(end-1)])
% % % title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% % % xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % % ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % % print('-dpng','-r300',['2D_Steady_Streamline_alpha_0']);
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% 
% % c_range = max(max(abs(omega_y)));
% plot(xp(:,end),zp(:,end),'.-k')
% plot(xTE(:,end),zTE(:,end),'.-g')
% if grd == 1
%     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% end
% % plot(xw,zw,'.-g')
% % plot(Xc,Zc,'xb')
% % streamline(Xf,Zf,U,W,Xf(:,1),Zf(:,1))
% 
% % bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
% % redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
% % vec = [bluevec; redvec(2:end,:)];
% % colormap(vec)
%     
% % pcolor(Xstar,Zstar,-omega_y)
% quiver(Xf,Zf,U,W,'k')
% shading interp
% % caxis(0.7*[-c_range c_range])
% axis equal
% axis([-c/2 + x_b(end-1) 3*c + x_b(end-1) -c/10 c])
% 
% axis([-c/2 + x_b(end-1) 3*c + x_b(end-1) -c + z_b(end-1) c + z_b(end-1)])
% % title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % print('-dpng','-r300',['2D_Wake_VF']);




% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xb')
% quiver(Xf,Yf,U,W,'b')
% axis equal
% axis([-c 4*c -c c])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xb')
% % quiver(Xf,Yf,u_d,w_d,'b')
% axis equal
% axis([-c 4*c -c c])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xb')
% quiver(Xf,Yf,u_s,w_s,'b')
% axis equal
% axis([-c 4*c -c c])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')


