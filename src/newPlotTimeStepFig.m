if i_t > 1 && PlotTimeStepFig == 1 
    FontSizeAx = 24;
    FontSizeLb = 32;
    afFigurePosition = [0 15 20 10];
    axespos = [0.2 0.2 0.8 0.8];
    ylabelpos = [-0.12 0.5];
    xlabelpos = [0.5 -0.20];

    % Create figure
    figure1 = fig('units','inches','width',12,'font','TimesNewRoman','fontsize',FontSizeAx);

    % Create axes
    axes1 = axes('Parent',figure1,'FontSize',FontSizeAx,'FontName','TimesNewRoman',...
        'PlotBoxAspectRatio',[3.25 1.5 1],...
        'Position',axespos,...
        'XTick',[-100:1:100],...
        'YTick',[-10:1:10],...
        'DataAspectRatio',[1 1 1]);

    hold(axes1,'on');

    % Create ylabel
    ylabel('$$z$$ (m)','Units','normalized','FontSize',FontSizeLb,...
        'FontName','TimesNewRoman',...
        'Rotation',90,...
        'Position',ylabelpos,...
        'Interpreter','latex');

    % Create xlabel
    xlabel('$$x$$ (m)','Units','normalized','FontSize',FontSizeLb,...
        'FontName','TimesNewRoman',...
        'Position',xlabelpos,...
        'Interpreter','latex');

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

    axis([-c/2 + x_b(2) 6*c + x_b(2) z_b(2) - 1.5*c z_b(2) + 1.5*c])

    % Create axes
    %axes2 = axes('Visible','on','Parent',figure1,...
    %    'PlotBoxAspectRatio',[2.5 1 4],...
    %    'Position',[0.67 0.37 0.3 0.8],...
    %    'XTick',[],...
    %    'YTick',[],...
    %    'DataAspectRatio',[1 1 1]);

    %box on
    %hold(axes2,'on');

    % Create plot
    %plot(xp,zp,'-k','linewidth',2)

    %axis([min(xp)-0.25 min(xp)+c+0.25 -0.375 0.375])
    if PrintFigures == 1
        figNum = [num2str(i_t,'%05d')];
        figName = sprintf('./movie/BEM2DFSI_%s.png',figNum);
        %figName = sprintf('./movie/BEM2DFSI_%s.eps',figNum);
        pause(1)
        print('-dpng','-r300',figName);
        %print('-depsc','-r300',figName);
        pause(0.25)
    else
        pause(1)
    end
    close(figure1)

    if PrintFigures == 1;
        figure2 = fig('units','inches','width',12,'font','TimesNewRoman','fontsize',FontSizeAx);

        % Create axes
        axes1 = axes('Visible','off','Parent',figure2,...
            'PlotBoxAspectRatio',[3.25 1.5 1],...
            'Position',axespos,...
            'XTick',[],...
            'YTick',[],...
            'DataAspectRatio',[1 1 1]);

        hold(axes1,'on');

        % Create plot
        plot(xp,zp,'-k','linewidth',2)

        axis([min(xp)-0.25 min(xp)+c+0.25 -0.375 0.375])

        figNum = [num2str(i_t,'%05d')];
        figName = sprintf('./movie/Body_%s.png',figNum);
        %figName = sprintf('./movie/Body_%s.eps',figNum);
        pause(1)
        print('-dpng','-r300',figName);
        %print('-depsc','-r300',figName);
        pause(0.25)
        close(figure2)
    end
end
   
