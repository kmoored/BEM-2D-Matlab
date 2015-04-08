if i_t > 1 && PlotTimeStepFig == 1 
    FontSizeAx = 24;
    FontSizeLb = 32;
    afFigurePosition = [0 15 20 10];
    axespos = [0.2 0.2 0.8 0.8];
    ylabelpos = [-0.12 0.5];
    xlabelpos = [0.5 -0.25];

    % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
    figHand = fig('units','inches','width',12,'font','TimesNewRoman','fontsize',FontSizeAx);
    %figHand = figure(1);
    subFigHand(1) = subplot(1,2,1);
    set(subFigHand(1),'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
    hold on
    axis equal
    plot(xp,zp,'-k','linewidth',2)

    if grd == 1
        plot(xp,-zp,'-k','linewidth',2)
        plot(xTE,-zTE,'.-b','linewidth',2)
        plot([min(xp)-10*c; min(xp)+10*c],[0 0],'-k','linewidth',2)
    end

    bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
    redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
    vec = [bluevec; redvec(2:end,:)];

    WakeCirc = GammaW'/max(abs(GammaW))/(1/1);
    WakeCirc(WakeCirc > 1) = 1;  
    WakeCirc(WakeCirc < -1) = -1;
    WakeCirc = WakeCirc + 1;
    WakeCirc = round(WakeCirc*100) + 1;
    WakeCirc(WakeCirc > 199) = 199;

    WakeCirc_w = -GammaW'/max(abs(GammaW))/(1/1);
    WakeCirc_w(WakeCirc_w > 1) = 1;  
    WakeCirc_w(WakeCirc_w < -1) = -1;
    WakeCirc_w = WakeCirc_w + 1;
    WakeCirc_w = round(WakeCirc_w*100) + 1;
    WakeCirc_w(WakeCirc_w > 199) = 199;

    WakeCirc_l = -muLump(1)'/max(abs(GammaW))/(1/1);
    WakeCirc_l(WakeCirc_l > 1) = 1;  
    WakeCirc_l(WakeCirc_l < -1) = -1;
    WakeCirc_l = WakeCirc_l + 1;
    WakeCirc_l = round(WakeCirc_l*100) + 1;
    WakeCirc_l(WakeCirc_l > 199) = 199;

    WakeCirc_lw = muLump(1)'/max(abs(GammaW))/(1/1);
    WakeCirc_lw(WakeCirc_lw > 1) = 1;  
    WakeCirc_lw(WakeCirc_lw < -1) = -1;
    WakeCirc_lw = WakeCirc_lw + 1;
    WakeCirc_lw = round(WakeCirc_lw*100) + 1;
    WakeCirc_lw(WakeCirc_lw > 199) = 199;


    for i_w = 1:wakeInd+1
        plot(xw(i_w),zw(i_w),'.-','color',vec(WakeCirc(i_w),:),'linewidth',2,'markersize',14)
    end
    if grd == 1 && i_t > 1
        for i_w = 1:wakeInd+1
            plot(xw(i_w),-zw(i_w),'.','color',vec(WakeCirc_w(i_w),:),'linewidth',2,'markersize',14)
        end
    end

    if i_t > Nlump*Nstep + 1
        plot(xl(1),zl(1),'.','color',vec(WakeCirc_l,:),'linewidth',2,'markersize',20)
        drawnow
        if grd == 1
            plot(xl(1),-zl(1),'.','color',vec(WakeCirc_lw,:),'linewidth',2,'markersize',20)

        end
    end

    %plot(xw,zw,'.k','linewidth',2,'markersize',14)
    %plot(xTE,zTE,'or')

    %hold off
    axis([-c/2 + x_b(2) 6*c + x_b(2) z_b(2) - 1.5*c z_b(2) + 1.5*c])

    set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
    %set(gca, 'Units', 'normalized', 'Position', axespos);
    % legend('d/c = \infty','d/c = 1/2','d/c = 1/4','Location','NorthWest')
    XTick = -100:1:100;
    YTick = -10:1:10;
    xlabel('$$x$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
    ylabel('$$z$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');

    %     print('-depsc','-r300',['Free_eps/2D_Pitch_Free_',num2str(i_t)]);
    %     print('-dpng','-r300',['Grd_png_free_lump4/2D_Pitch_Grd_',num2str(i_t)]);
    %     pause(0.01);
    
    subFigHand(2) = subplot(1,2,2);
    plot(xp,zp,'-k','linewidth',2)
    axis equal
    axis([min(xp)-0.125 min(xp)+c+0.125 -0.25 0.25])
    axis off
    %box off
    
    %hold off
    figNum = [num2str(i_t,'%05d')];
    %figName = sprintf('./movie/BEM2DFSI_%s.png',figNum);
    figName = sprintf('./movie/BEM2DFSI_%s.eps',figNum);
    pause(0.25)
    %print('-dpng','-r300',figName);
    print('-depsc',figName);
    pause(0.25)
    %close(figHand)
end
   
