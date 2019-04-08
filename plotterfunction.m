%% This script just plots the result for G (the current distribution) 
%
% 1. Compute the magnetic field
% 2. plot the currents in 2D
% 3. plot the currents and the magnetic field in z-direction in 3D
% 4. plot a cut trough the current distribution along the y-direction
% 5. plot energy history

function plotterfunction(G,geometrymask,geometry,E_vec,whole_current,Ha)


name = 'hole';

my = 4*pi*1e-07;

dX = geometry.dX;
dY = geometry.dY;
deltaX = geometry.deltaX;
deltaY = geometry.deltaY;
gridpointX = geometry.gridpointX;
gridpointY = geometry.gridpointY;
space = geometry.space;
nospace = geometry.nospace;
W = 10^(-5);

[Bz, Kx, Ky, K_abs_array] = calculator(G,geometrymask,geometry);

if  whole_current~=0
    K_abs_array = K_abs_array*W/whole_current;
    Bz = W*Bz/(my*whole_current);
    curr_label = 'K/(I_{a}/W)';
    field_label = 'B_{z}/(\mu_{0}I_{a}/W)';
else    
    K_abs_array = K_abs_array/Ha; % 10^5 korrektur ... ?
    Bz = Bz/(my*Ha); % 10^5 korrektur ... ?
    curr_label = 'K/(H_{a})';
    field_label = 'B_{z}/(\mu_{0}H_{a})';
end


%% This part here plots the currents an their densities in a cell (2D-plot)

measureX_1 = (deltaX/2:deltaX:dX-deltaX/2);
measureY_1 = (deltaY/2:deltaY:dY-deltaY/2);
measureY_1 = fliplr(measureY_1);

[X1,Y1] = meshgrid(measureX_1,measureY_1);

curr2d = figure;
    curr2d.Position = [100   514   570   450];
%     set(gcf,'color','white')
    imagesc([min(measureX_1) max(measureX_1)], [min(measureY_1) max(measureY_1)], flipud(K_abs_array)); 
    hold on;
    
    c = hot(1000); c = flipud(c(250:900,:));
    c = [ ones(1,3); c];
    colormap(c);
    c = colorbar;
    c.Label.String = curr_label;
    c.Label.Rotation = 0;
    c.Label.Units = 'normalized';
    c.Label.Position = [-0.1 1.075 0];
    c.Label.FontSize = 12;
    
    axis([0 dX 0 dY]);
    ax = gca;
    ax.YTick = [ 0.5 1 1.5 ]*1E-5;
    ax.XTick = [ 0.5 1 1.5 ]*1E-5;
    xlabel('x [m]');
    ylabel('y [m]');   
    set(gca,'ydir','normal');
    hold off;
    stream = streamslice(X1,Y1,Kx,-1*Ky,1);
    set( stream, 'Color', [0.4 0.4 0.4] );
    title ('Current distribution K_{abs}');

% export_fig(curr2d, [name, '_current_currcut.pdf'], '-pdf');
    
    
%% Here we plot the magnetic field over the currents (3D-plot)

measureX_2 = (deltaX:deltaX:dX-deltaX);
measureY_2 = (deltaY:deltaY:dY-deltaY);
measureY_2 = fliplr(measureY_2);
[X2,Y2] = meshgrid(measureX_2,measureY_2);

field3d = figure;
    field3d.Position = [800   514   700   450];
%     set(gcf,'color','white')
    surf(X2,Y2,Bz);
    view(45,45);
    c = colorbar;
    
    c.Label.String = field_label;
    c.Label.Rotation = 0;
    c.Label.Units = 'normalized';
    c.Label.Position = [-0.1 1.075 0];
    c.Label.FontSize = 12;
    
    axis([0 dX 0 dY]);
    ax = gca;
    ax.YTick = [ 0.5 1 1.5 ]*1E-5;
    ax.XTick = [ 0.5 1 1.5 ]*1E-5;
    xlabel('x [m]');
    ylabel('y [m]');
    title ('Magnetic field at z=0');

% export_fig(field3d, [name, '_current_fieldcut.pdf'], '-pdf');


end

function [Bz, Kx, Ky, K_abs_array] = calculator(G,geometrymask,geometry)

my = 4*pi*1e-07; %%a physical constant ...

dX = geometry.dX;
dY = geometry.dY;
deltaX = geometry.deltaX;
deltaY = geometry.deltaY;
gridpointX = geometry.gridpointX;
gridpointY = geometry.gridpointY;
space = geometry.space;
nospace = geometry.nospace;

%% Some computations

%%Computing the currents
Kx = zeros(gridpointY-1,gridpointX-1);       
Ky = zeros(gridpointY-1,gridpointX-1);
K_abs_array = zeros(gridpointY-1,gridpointX-1);
for i = 1:1:gridpointY-1
for j = 1:1:gridpointX-1
    if (geometrymask(i,j) ~= nospace) && (geometrymask(i+1,j) ~= nospace) ...
       && (geometrymask(i,j+1) ~= nospace) && (geometrymask(i+1,j+1) ~= nospace)
        %%(Kx,Ky) = (dyG,-dxG) because of K = nabla X (g(x,y)*ez)
        Kx(i,j) = (0.5/deltaY)*((G(i+1,j)-G(i,j))+(G(i+1,j+1)-G(i,j+1)));
        Ky(i,j) = (-1)*(0.5/deltaX)*((G(i,j+1)-G(i,j))+(G(i+1,j+1)-G(i+1,j)));
        K_abs_array(i,j) = sqrt(Kx(i,j)*Kx(i,j)+Ky(i,j)*Ky(i,j));
    end
end
end

%%Vector potential (A = (Ax,Ay,Az))
Ax = zeros(gridpointY-1,gridpointX-1);
Ay = zeros(gridpointY-1,gridpointX-1);
x = deltaX/2:deltaX:dX-deltaX/2;
y = deltaY/2:deltaY:dY-deltaY/2;
x2 = x;
y2 = y;
for i = 1:1:gridpointY-1
for j = 1:1:gridpointX-1
    sumAx = 0; sumAy = 0;
    for i2 = 1:1:gridpointY-1
    for j2 = 1:1:gridpointX-1

        dummy_Ax = (1/(4*pi))*deltaX*deltaY*geometry.thickness*my*Kx(i2,j2)/(abs((x(j)-x2(j2))^2+(y(i)-y2(i2))^2));
        dummy_Ay = (1/(4*pi))*deltaX*deltaY*geometry.thickness*my*Ky(i2,j2)/(abs((x(j)-x2(j2))^2+(y(i)-y2(i2))^2));
        if isinf(dummy_Ax) || isnan(dummy_Ax)
            dummy_Ax = 0;
        end
        if isinf(dummy_Ay) || isnan(dummy_Ay)
            dummy_Ay = 0;
        end
        sumAx = sumAx + dummy_Ax;
        sumAy = sumAy + dummy_Ay;
    end
    end
    Ax(i,j) = sumAx;
    Ay(i,j) = sumAy;
end
end

%%And the magnetic field
Bz = zeros(gridpointY-2,gridpointX-2);
for i = 1:1:gridpointY-2
for j = 1:1:gridpointX-2
    %%(Bx,By,Bz) = (,0,dxAy-dyAx)
    Bz(i,j) = (0.5/deltaX)*((Ay(i,j+1)-Ay(i,j))+(Ay(i+1,j+1)-Ay(i+1,j)))-...
              (0.5/deltaY)*((Ax(i+1,j)-Ax(i,j))+(Ax(i+1,j+1)-Ax(i,j+1)));
end
end


end
   