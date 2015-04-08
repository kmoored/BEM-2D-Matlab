%% Free-swimming calculations

F_t(i_t) = Fx;   
if (free == 1 && i_t > 1)
    a_b = Fx/M;       
    Q0(1,i_t+1) = a_b*delT + Q0(1,i_t);

    x_b(1) = (Q0(1,i_t + 1) + Q0(1,i_t))*delT/2 + x_b(2); 


    if PlotVelFig == 1 && (rem(i_t,5) == 0)
        if i_t > 10
            close(velh(i_t-10))
        end

        FontSizeAx = 24;
        FontSizeLb = 32;
        afFigurePosition = [0 0 20 10];
        axespos = [0.15 0.15 0.84 0.75];
        ylabelpos = [-0.12 0.5];
        xlabelpos = [0.5 -0.25];

        velh(i_t) = figure;
        set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
        set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

        hold on
            plot([0 Ncyc/f],U_swim*[1 1],'-k','linewidth',1)
            plot(N_a/f*[1 1],[0 1.5*max(-Q0(1,1:i_t))],'--k','linewidth',1)
            plot(((1:i_t) - 1)*delT,-Q0(1,1:i_t)','-b','linewidth',2)
        hold off

        set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
        set(gca, 'Units', 'normalized', 'Position', axespos);
        xlabel('$$t/T$$, sec','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
        ylabel('$$U(t)$$','Interpreter', 'LaTeX','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos);
        axis([0 Ncyc/f 0 1.5*max(-Q0(1,1:i_t))])
    end
%         subplot(2,1,2)
%         plot((1:Ncyc*Nstep+1)*delT,D_visc,'-r','linewidth',2)
%         xlabel('$$t$$, sec','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
%         ylabel('Viscous Drag','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
%         axis([0 Ncyc/f 0 1.5*max(D_visc)])
else
    x_b(1) = Q0(1,i_t)*delT + x_b(2);
end
z_b(1) = Q0(2,i_t)*delT + z_b(2);


%     close(wbar);
% delTime(i_t) = toc;
% EstT = delTime(i_t)*(Ncyc*Nstep+1 - i_t);
% EstThour = floor(EstT/3600);
% EstTmin = floor((EstT - EstThour*3600)/60);
% EstTsec = round(EstT - EstThour*3600 - EstTmin*60);
%     sprintf(['Calculating...Time Remaining: ',num2str(EstThour),' hrs ',num2str(EstTmin),' mins ',num2str(EstTsec),' secs']);

%     wbar = waitbar(i_t/(Ncyc*Nstep+1),['Calculating...Time Remaining: ',num2str(EstThour),' hrs ',num2str(EstTmin),' mins ',num2str(EstTsec),' secs'],'Position',[(scrsz(3) - scrsz(4)/3) 0 scrsz(4)/3 1/16*scrsz(4)]);

if SaveData == 1           
    Data = [wakeInd,Q0(1,i_t),Q0(2,i_t),x_b(2),z_b(2),a_b,Fx,Fz,Pow,L,T,D_visc,Gamma,Cf,Cl,Ct,Cpow,...
        Fx_s,Fx_us,Fz_s,Fz_us,Pow_s,Pow_us,L_s,L_us,T_s,T_us,Cl_s,Cl_us,...
        Ct_s,Ct_us,Cpow_s,Cpow_us];
    fprintf(fid_Data,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',Data');

    PanelProp = [xp_0(1:end-1)'; zp_0(1:end-1)';xp(1:end-1)';zp(1:end-1)';...
        Vc(:,1)';Vc(:,2)';xc';zc';vt(:,1)';vt(:,2)';vn(:,1)';vn(:,2)';...
        dL';Xc';Zc';sigma';mu(:,1)';Qp';Qt';Cp_s';Cp_us';Cp';dFshear']';
    fprintf(fid_PanelProp,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',PanelProp');

    WakeProp = [xw xl(1); zw zl(1);vtTE' vtw' vtlw(1,:)';...
        vnTE' vnw' vnlw(1,:)'; muTE(1) muW muLump(1); 
        GammaW -muLump(1)]';
    fprintf(fid_WakeProp,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',WakeProp');
end

%fprintf([num2str(round((i_t/(Ncyc*Nstep+1))*100)),'%% complete\n'],'FontName','TimesNewRoman')
