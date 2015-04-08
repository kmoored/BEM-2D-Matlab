function PlotPanelsMatrix(t,pc,pcI,Vc,pcTE,pcTEI,pcw,pcwI,C1,cpts,cptsI,cptsTE,cptsTEI,cptsw,cptswI,vc,vcI,vn,vnI,vt,vtI,vcTE,vcTEI,vnTE,vnTEI,vtTE,vtTEI,vcw,vcwI,vnw,vnwI,vtw,vtwI,Npanels,Nspanels,Cp,tau,SurfaceColor)

[~,Nw] = size(pcw);

% Decomposing panel corner points into x, y and z coordinates.
pc_x = pc(1:3:12,:);
pc_y = pc(2:3:12,:);
pc_z = pc(3:3:12,:);

pcI_x = pcI(1:3:12,:);
pcI_y = pcI(2:3:12,:);
pcI_z = pcI(3:3:12,:);

pcTE_x = pcTE(1:3:12,:);
pcTE_y = pcTE(2:3:12,:);
pcTE_z = pcTE(3:3:12,:);

pcTEI_x = pcTEI(1:3:12,:);
pcTEI_y = pcTEI(2:3:12,:);
pcTEI_z = pcTEI(3:3:12,:);

if t > 0
    pcw_x = pcw(1:3:12,:);
    pcw_y = pcw(2:3:12,:);
    pcw_z = pcw(3:3:12,:);
    
    pcwI_x = pcwI(1:3:12,:);
    pcwI_y = pcwI(2:3:12,:);
    pcwI_z = pcwI(3:3:12,:);
end

span = 2*max(max(abs(pc_y)));

% Cycle thru the panels
for k = 1:2*Npanels
    
    % Cycle thru the panel edges.
    for i = 1:4
        
        % Use the connectivity to find coordinate pairs.
        
        % Body panels.
        eleX = [(C1(i,:) == 1)*pc_x(:,k);(C1(i,:) == -1)*pc_x(:,k)];
        eleY = [(C1(i,:) == 1)*pc_y(:,k);(C1(i,:) == -1)*pc_y(:,k)];
        eleZ = [(C1(i,:) == 1)*pc_z(:,k);(C1(i,:) == -1)*pc_z(:,k)];
%         plot3(eleX,eleY,eleZ,'-b', 'LineWidth',1);
        
        % Body image panels.
        eleX = [(C1(i,:) == 1)*pcI_x(:,k);(C1(i,:) == -1)*pcI_x(:,k)];
        eleY = [(C1(i,:) == 1)*pcI_y(:,k);(C1(i,:) == -1)*pcI_y(:,k)];
        eleZ = [(C1(i,:) == 1)*pcI_z(:,k);(C1(i,:) == -1)*pcI_z(:,k)];
%         plot3(eleX,eleY,eleZ,'-b', 'LineWidth',1);
        
    end
    % Plotting collocation points.
%     plot3(cpts(1,k),cpts(2,k),cpts(3,k),'xr')
%     plot3(cptsI(1,k),cptsI(2,k),cptsI(3,k),'xr')
    
    % Plotting directional vectors.
    val = 10/span;
%     quiver3(cpts(1,k),cpts(2,k),cpts(3,k),1/val*vc(1,k),1/val*vc(2,k),1/val*vc(3,k),'k')
%     quiver3(cpts(1,k),cpts(2,k),cpts(3,k),1/val*vn(1,k),1/val*vn(2,k),1/val*vn(3,k),'r')
%     quiver3(cpts(1,k),cpts(2,k),cpts(3,k),1/val*vt(1,k),1/val*vt(2,k),1/val*vt(3,k),'k')
%     quiver3(cpts(1,k),cpts(2,k),cpts(3,k),1/val*Vc(1,k),1/val*Vc(2,k),1/val*Vc(3,k),'b')
%     
%     quiver3(cptsI(1,k),cptsI(2,k),cptsI(3,k),1/val*vcI(1,k),1/val*vcI(2,k),1/val*vcI(3,k),'k')
%     quiver3(cptsI(1,k),cptsI(2,k),cptsI(3,k),1/val*vnI(1,k),1/val*vnI(2,k),1/val*vnI(3,k),'k')
%     quiver3(cptsI(1,k),cptsI(2,k),cptsI(3,k),1/val*vtI(1,k),1/val*vtI(2,k),1/val*vtI(3,k),'k')
   
end

for k = 1:Nspanels

    % Cycle thru the panel edges.
    for i = 1:4

        % Use the connectivity to find coordinate pairs.

        % TE panels.
        eleX = [(C1(i,:) == 1)*pcTE_x(:,k);(C1(i,:) == -1)*pcTE_x(:,k)];
        eleY = [(C1(i,:) == 1)*pcTE_y(:,k);(C1(i,:) == -1)*pcTE_y(:,k)];
        eleZ = [(C1(i,:) == 1)*pcTE_z(:,k);(C1(i,:) == -1)*pcTE_z(:,k)];
        plot3(eleX,eleY,eleZ,'-k', 'LineWidth',1);
        
        % TE image panels.
        eleX = [(C1(i,:) == 1)*pcTEI_x(:,k);(C1(i,:) == -1)*pcTEI_x(:,k)];
        eleY = [(C1(i,:) == 1)*pcTEI_y(:,k);(C1(i,:) == -1)*pcTEI_y(:,k)];
        eleZ = [(C1(i,:) == 1)*pcTEI_z(:,k);(C1(i,:) == -1)*pcTEI_z(:,k)];
        plot3(eleX,eleY,eleZ,'-k', 'LineWidth',1);

    end

    % Plotting collocation points.
%     plot3(cptsTE(1,k),cptsTE(2,k),cptsTE(3,k),'xr')
%     plot3(cptsTEI(1,k),cptsTEI(2,k),cptsTEI(3,k),'xr')

    % Plotting directional vectors.
%         val = 10/span;
%         quiver3(cptsTE(1,k),cptsTE(2,k),cptsTE(3,k),1/val*vnTE(1,k),1/val*vnTE(2,k),1/val*vnTE(3,k),'k')
%         quiver3(cptsTE(1,k),cptsTE(2,k),cptsTE(3,k),1/val*vcTE(1,k),1/val*vcTE(2,k),1/val*vcTE(3,k),'k')
%         quiver3(cptsTE(1,k),cptsTE(2,k),cptsTE(3,k),1/val*vtTE(1,k),1/val*vtTE(2,k),1/val*vtTE(3,k),'k')
%         
%         quiver3(cptsTEI(1,k),cptsTEI(2,k),cptsTEI(3,k),1/val*vnTEI(1,k),1/val*vnTEI(2,k),1/val*vnTEI(3,k),'k')
%         quiver3(cptsTEI(1,k),cptsTEI(2,k),cptsTEI(3,k),1/val*vcTEI(1,k),1/val*vcTEI(2,k),1/val*vcTEI(3,k),'k')
%         quiver3(cptsTEI(1,k),cptsTEI(2,k),cptsTEI(3,k),1/val*vtTEI(1,k),1/val*vtTEI(2,k),1/val*vtTEI(3,k),'k')
end

if t > 0
    for k = 1:Nw

        % Cycle thru the panel edges.
        for i = 1:4

            % Use the connectivity to find coordinate pairs.

            % Wake panels.
            eleX = [(C1(i,:) == 1)*pcw_x(:,k);(C1(i,:) == -1)*pcw_x(:,k)];
            eleY = [(C1(i,:) == 1)*pcw_y(:,k);(C1(i,:) == -1)*pcw_y(:,k)];
            eleZ = [(C1(i,:) == 1)*pcw_z(:,k);(C1(i,:) == -1)*pcw_z(:,k)];
            plot3(eleX,eleY,eleZ,'-k', 'LineWidth',1);

            % Wake image panels.
            eleX = [(C1(i,:) == 1)*pcwI_x(:,k);(C1(i,:) == -1)*pcwI_x(:,k)];
            eleY = [(C1(i,:) == 1)*pcwI_y(:,k);(C1(i,:) == -1)*pcwI_y(:,k)];
            eleZ = [(C1(i,:) == 1)*pcwI_z(:,k);(C1(i,:) == -1)*pcwI_z(:,k)];
            plot3(eleX,eleY,eleZ,'-k', 'LineWidth',1);

        end

        % Plotting collocation points.
%             plot3(cptsw(1,k,kw),cptsw(2,k,kw),cptsw(3,k,kw),'xr')
%             plot3(cptswI(1,k,kw),cptswI(2,k,kw),cptswI(3,k,kw),'xr')

        % Plotting directional vectors.
%             val = 10/span;
%             quiver3(cptsw(1,k,kw),cptsw(2,k,kw),cptsw(3,k,kw),1/val*vnw(1,k,kw),1/val*vnw(2,k,kw),1/val*vnw(3,k,kw),'k')
%             quiver3(cptsw(1,k,kw),cptsw(2,k,kw),cptsw(3,k,kw),1/val*vcw(1,k,kw),1/val*vcw(2,k,kw),1/val*vcw(3,k,kw),'k')
%             quiver3(cptsw(1,k,kw),cptsw(2,k,kw),cptsw(3,k,kw),1/val*vtw(1,k,kw),1/val*vtw(2,k,kw),1/val*vtw(3,k,kw),'k')
%             
%             quiver3(cptswI(1,k,kw),cptswI(2,k,kw),cptswI(3,k,kw),1/val*vnwI(1,k,kw),1/val*vnwI(2,k,kw),1/val*vnwI(3,k,kw),'k')
%             quiver3(cptswI(1,k,kw),cptswI(2,k,kw),cptswI(3,k,kw),1/val*vcwI(1,k,kw),1/val*vcwI(2,k,kw),1/val*vcwI(3,k,kw),'k')
%             quiver3(cptswI(1,k,kw),cptswI(2,k,kw),cptswI(3,k,kw),1/val*vtwI(1,k,kw),1/val*vtwI(2,k,kw),1/val*vtwI(3,k,kw),'k')
    end
end

if SurfaceColor == 1
    cdata = -Cp'; 
    colormap(jet)
    caxis([-1 2]);
else
    cdata = abs(tau)'; 
    colormap(jet)
    caxis([0 max(abs(tau))]);
end

ptch = patch(pc_x,pc_y,pc_z,'c');
ptchI = patch(pcI_x,pcI_y,pcI_z,'c');
set(ptch,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled','EdgeColor','none')
set(ptchI,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled','EdgeColor','none')

if t > 0
    ptchw = patch(pcw_x,pcw_y,pcw_z,'g');
    ptchwI = patch(pcwI_x,pcwI_y,pcwI_z,'g');
    set(ptchw,'FaceColor','c','FaceAlpha',0.25,'FaceLighting','Phong')
    set(ptchwI,'FaceColor','c','FaceAlpha',0.25,'FaceLighting','Phong')
end
% colorbar
axis equal


