function PlotGrd(Ncyc,Nstep,delT,Q0,c_b,b,c_r,C1)

Ny = 50;
Nx = 50;

xminus = (Ncyc*Nstep + 2)*delT*Q0(1) - (c_b*b - c_r)/2*1.05 - 4*b;
xplus = -(Ncyc*Nstep + 2)*delT*Q0(1) + (c_b*b - c_r)/2*1.05 + 4*b;
yminus = -6*b;
yplus = 6*b;
x = linspace(xminus,xplus,Nx);
y = linspace(yminus,yplus,Ny);

[X,Y] = meshgrid(x,y);
Z = zeros(Nx,Ny);

k = 1;
for j = 1:Ny - 1
    for i = 1:Nx - 1
        pc(:,k) = [X(i,j);Y(i,j);Z(i,j);  X(i,j+1);Y(i,j+1);Z(i,j+1);  X(i+1,j+1);Y(i+1,j+1);Z(i+1,j+1);  X(i+1,j);Y(i+1,j);Z(i+1,j)];
        
        k = k + 1;
    end
end
GrdPanels = k-1;

[pc_x,pc_y,pc_z] = Decompose(pc);

% Cycle thru the panels
for k = 1:GrdPanels
    
    % Cycle thru the panel edges.
    for i = 1:4
        
        % Use the connectivity to find coordinate pairs.
        eleX = [(C1(i,:) == 1)*pc_x(:,k);(C1(i,:) == -1)*pc_x(:,k)];
        eleY = [(C1(i,:) == 1)*pc_y(:,k);(C1(i,:) == -1)*pc_y(:,k)];
        eleZ = [(C1(i,:) == 1)*pc_z(:,k);(C1(i,:) == -1)*pc_z(:,k)];
        plot3(eleX,eleY,eleZ,'-k', 'LineWidth',1);
    end
end

