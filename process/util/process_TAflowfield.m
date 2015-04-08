function process_TAflowfield(folder,savefilename)

%% Load Processed Data
load([folder,'/Processed',savefilename,'.mat']);
load([folder,'/Flowfield',savefilename,'.mat']); 

%% Plotting

for i_t = Ncyc*Nstep+1
    figure;
    FontSizeAx = 24;
    afFigurePosition = [15 7 25 15];

    set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
    set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

    hold on
    axis equal


    bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
    redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
    vec = [bluevec; redvec(2:end,:)];
    colormap(vec)

%     Vmag = mean(sqrt(Ut.^2 + Wt.^2),3);
    TA_W = mean(Wt(:,:,:),3)
    TA_U = mean(Ut(:,:,:),3)

    pcolor(Xf,Zf,TA_W)
%     quiver(Xf,Zf,mean(Ut,3),mean(Wt,3),'k')
    % quiver(Xf,Zf,Uinf*ones(nz,nx),mean(Wt(:,:,1:151)-Winf,3),'k')
    plot(xp_0 + x_b(end) ,zp_0 + d_c*c,'-k','linewidth',2)
    hold off

    %         if grd == 1
    %             plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
    %         end

    shading interp
    colormap jet
    caxis(0.3*[-1 1])
    colorbar
    axis([-c/4 + x_b(i_t) 6*c + x_b(i_t) -c c])
    print('-dpng','-r600',[folder,'/TA_W_',savefilename,'.png']);
    
    figure
    FontSizeAx = 24;
    afFigurePosition = [15 7 25 15];

    set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
    set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

    hold on
    axis equal
    
     pcolor(Xf,Zf,TA_U)
%     quiver(Xf,Zf,mean(Ut,3),mean(Wt,3),'k')
    % quiver(Xf,Zf,Uinf*ones(nz,nx),mean(Wt(:,:,1:151)-Winf,3),'k')
    plot(xp_0 + x_b(end) ,zp_0 + d_c*c,'-k','linewidth',2)
    hold off

    %         if grd == 1
    %             plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
    %         end

    shading interp
    colormap jet
    caxis(1.5*[1/3 1])
    colorbar
    axis([-c/4 + x_b(i_t) 6*c + x_b(i_t) -c c])
    print('-dpng','-r600',[folder,'/TA_U_',savefilename,'.png']);
    %         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
    %         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
    %         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

    %         VortMov(i_t-1) = getframe;
    %         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);
            

    %             print('-dpng','-r300',['GrdData/Vorticity/GrdEffect_St',num2str(St*1000),'_d_c',num2str(d_c*100),'_A_c',num2str(A_c),'_',num2str(figInd),'.png']);
end