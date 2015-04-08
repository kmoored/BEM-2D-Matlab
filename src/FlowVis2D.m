clear
clc

% d_c_vec = [100 1 5/6 4/6 3/6]';
% St_vec = [0.025:0.025:0.5]';
St = 0.3;
d_c = 5/6;
A_c = 0.25;

%% Loading data
load(['Data_GrdEffect2D/Pitch_St',num2str(St*100),...
    '_A_c',num2str(A_c*100),...
    '_d_c',num2str(d_c*100),'.mat']);

%% Airfoil plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [0 2 30 20];
axespos = [0.15 0.2 0.84 0.8];
ylabelpos = [-0.12 0.5];
xlabelpos = [0.5 -0.25];


% Plotting airfoil with LE vortex sheets and the TE vortex sheet.
fighand = figure;
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
axis equal

plot(Xc(:,end),Zc(:,end),'xr','linewidth',2)
plot(Xc(1,end),Zc(1,end),'xb','linewidth',2)
plot(Xc(end,end),Zc(end,end),'xg','linewidth',2)
plot(xp(:,end),zp(:,end),'.-k','linewidth',2)
plot(xTE(:,end),zTE(:,end),'.-b','linewidth',2)
% quiver(xc(:,end),zc(:,end),vn(:,1,end),vn(:,2,end),'b','linewidth',2)
val = 1/10;

if grd == 1
%         plot(xp_2(:,i_t),zp_2(:,i_t),'-k','linewidth',2)
%         plot(xTE_2(:,i_t),zTE_2(:,i_t),'.-b','linewidth',2)
    plot([min(xp(:,end))-10*c; min(xp(:,end))+10*c],[0 0],'-k','linewidth',4)
%         if t > 0
%             plot(xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),'.-b','linewidth',2,'markersize',14)
%         end
end

if t > 0
    if LES == 1
        plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
        plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)

        for j = i_t:-1:2
            if muLEt(Ncyc*Nstep + 2 - j) == 0
            else
                plot(xpt_LES(:,Ncyc*Nstep + 2 - j),zpt_LES(:,Ncyc*Nstep + 2 - j),'.-g','linewidth',2,'markersize',14)
%                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vtLEt(Ncyc*Nstep + 1 - j,1),val*vtLEt(Ncyc*Nstep + 1 - j,2),'b')
%                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vnLEt(Ncyc*Nstep + 1 - j,1),val*vnLEt(Ncyc*Nstep + 1 - j,2),'k')

            end
        end

        for j = i_t:-1:2
            if muLEb(Ncyc*Nstep + 2 - j) == 0  
            else
                plot(xpb_LES(:,Ncyc*Nstep + 2 - j),zpb_LES(:,Ncyc*Nstep + 2 - j),'.-r','linewidth',2,'markersize',14)
%                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vtLEb(Ncyc*Nstep + 1 - j,1),val*vtLEb(Ncyc*Nstep + 1 - j,2),'b')
%                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vnLEb(Ncyc*Nstep + 1 - j,1),val*vnLEb(Ncyc*Nstep + 1 - j,2),'k')

            end
        end
    end
%         plot(xp(Stagpt,i_t),zp(Stagpt,i_t),'xk','linewidth',2)

    bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
    redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
    vec = [bluevec; redvec(2:end,:)];
    vec = vec(end:-1:1,:);

    WakeCirc = GammaW'/max(GammaW)/(1/2);
    WakeCirc(WakeCirc > 1) = 1;  
    WakeCirc(WakeCirc < -1) = -1;
    WakeCirc = WakeCirc + 1;
    WakeCirc = round(WakeCirc*100) + 1;
    WakeCirc(WakeCirc > 199) = 199;

    for i_w = 1:Ncyc*Nstep+1
        plot(xw(i_w),zw(i_w),'.','color',vec(WakeCirc(i_w),:),'linewidth',2,'markersize',14)
    end
end

axis equal
%     axis([-c/2 + x_b(i_t) 2*c + x_b(i_t) -1/2*c + z_b(i_t) 1/2*c + z_b(i_t)])
axis([-c/2 + x_b(end-1) 7*c + x_b(end-1) z_b(end-1) - 1*c z_b(end-1) + 1*c])

set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$x$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$z$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');

