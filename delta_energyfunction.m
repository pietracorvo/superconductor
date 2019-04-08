%%This script computes the new energy of a cell

function [dE] = delta_energyfunction(G,dG,Ha,Lambda,Nconst,gradGY_vec,gradGX_vec,neighbourmask,geometry)

%%This functon computes the energies just for 3 x 3 grid fragments of G,
%%so just for a certain gridpoint (i,j) that has been variated 
%%with the 8 adjacent points of G, which stay the same.

my = 4*pi*1e-07;                   %%a physical constant!           

deltaX = geometry.deltaX;
deltaY = geometry.deltaY;
real_cellnumber = geometry.real_cellnumber;
nospace = geometry.nospace;


%% Variate the G

G_old = G;
G(2,2) = G(2,2) + dG;

%%save old_ data
old_gradGY_vec = gradGY_vec;
old_gradGX_vec = gradGX_vec;
                
                
%% Compute the new gradients

old_gradGY = zeros(2,2);
old_gradGX = zeros(2,2);
new_gradGY = zeros(2,2);
new_gradGX = zeros(2,2);

if neighbourmask(1,1) ~= nospace
    k = neighbourmask(1,1);
    old_gradGY(1,1) = (0.5/deltaY) * (G_old(2,1)-G_old(1,1)+G_old(2,2)-G_old(1,2));
    old_gradGX(1,1) = (0.5/deltaX) * (G_old(1,2)-G_old(1,1)+G_old(2,2)-G_old(2,1));
    new_gradGY(1,1) = (0.5/deltaY) * (G(2,1)-G(1,1)+G(2,2)-G(1,2));
    new_gradGX(1,1) = (0.5/deltaX) * (G(1,2)-G(1,1)+G(2,2)-G(2,1));
    gradGY_vec(k) = new_gradGY(1,1);
    gradGX_vec(k) = new_gradGX(1,1);
end
if neighbourmask(1,2) ~= nospace
    k = neighbourmask(1,2);
    old_gradGY(1,2) = (0.5/deltaY) * (G_old(2,2)-G_old(1,2)+G_old(2,3)-G_old(1,3));
    old_gradGX(1,2) = (0.5/deltaX) * (G_old(1,3)-G_old(1,2)+G_old(2,3)-G_old(2,2));
    new_gradGY(1,2) = (0.5/deltaY) * (G(2,2)-G(1,2)+G(2,3)-G(1,3));
    new_gradGX(1,2) = (0.5/deltaX) * (G(1,3)-G(1,2)+G(2,3)-G(2,2));
    gradGY_vec(k) = new_gradGY(1,2);
    gradGX_vec(k) = new_gradGX(1,2);
end
if neighbourmask(2,1) ~= nospace
    k = neighbourmask(2,1);
    old_gradGY(2,1) = (0.5/deltaY) * (G_old(3,1)-G_old(2,1)+G_old(3,2)-G_old(2,2));
    old_gradGX(2,1) = (0.5/deltaX) * (G_old(2,2)-G_old(2,1)+G_old(3,2)-G_old(3,1));
    new_gradGY(2,1) = (0.5/deltaY) * (G(3,1)-G(2,1)+G(3,2)-G(2,2));
    new_gradGX(2,1) = (0.5/deltaX) * (G(2,2)-G(2,1)+G(3,2)-G(3,1));
    gradGY_vec(k) = new_gradGY(2,1);
    gradGX_vec(k) = new_gradGX(2,1);
end
if neighbourmask(2,2) ~= nospace
    k = neighbourmask(2,2);
    old_gradGY(2,2) = (0.5/deltaY) * (G_old(3,2)-G_old(2,2)+G_old(3,3)-G_old(2,3));
    old_gradGX(2,2) = (0.5/deltaX) * (G_old(2,3)-G_old(2,2)+G_old(3,3)-G_old(3,2));
    new_gradGY(2,2) = (0.5/deltaY) * (G(3,2)-G(2,2)+G(3,3)-G(2,3));
    new_gradGX(2,2) = (0.5/deltaX) * (G(2,3)-G(2,2)+G(3,3)-G(3,2));
    gradGY_vec(k) = new_gradGY(2,2);
    gradGX_vec(k) = new_gradGX(2,2);
end    


%% Compute new parts of Eint

old_Eint_array = zeros(real_cellnumber,8);
new_Eint_array = zeros(real_cellnumber,8);

if neighbourmask(1,1) ~= nospace
    k = neighbourmask(1,1);
    old_Eint_array(:,1) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,1);
    old_Eint_array(:,2) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,1); 
    new_Eint_array(:,1) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,1);
    new_Eint_array(:,2) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,1);
end

if neighbourmask(1,2) ~= nospace
    k = neighbourmask(1,2);
    old_Eint_array(:,3) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,2);
    old_Eint_array(:,4) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,2); 
    new_Eint_array(:,3) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,2);
    new_Eint_array(:,4) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,2);
end

if neighbourmask(2,1) ~= nospace
    k = neighbourmask(2,1);
    old_Eint_array(:,5) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,3);
    old_Eint_array(:,6) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,3); 
    new_Eint_array(:,5) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,3);
    new_Eint_array(:,6) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,3);
end

if neighbourmask(2,2) ~= nospace
    k = neighbourmask(2,2);
    old_Eint_array(:,7) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,4);
    old_Eint_array(:,8) = 0.5 * my * (old_gradGY_vec*old_gradGY_vec(k)+old_gradGX_vec*old_gradGX_vec(k)).*Nconst(:,4); 
    new_Eint_array(:,7) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,4);
    new_Eint_array(:,8) = 0.5 * my * (gradGY_vec*gradGY_vec(k)+gradGX_vec*gradGX_vec(k)).*Nconst(:,4);
end

%%some values are counted twice, we correct that here ...
%%That corrects multiple counting of values in Eint_array 
if neighbourmask(1,1) ~= nospace  %%These 4 if's are for double counting in row and column for one k ...
    k = neighbourmask(1,1);
    new_Eint_array(k,2) = 0;
    old_Eint_array(k,2) = 0;
end
if neighbourmask(1,2) ~= nospace
    k = neighbourmask(1,2);
    old_Eint_array(k,4) = 0;
    new_Eint_array(k,4) = 0;
end
if neighbourmask(2,1) ~= nospace
    k = neighbourmask(2,1);
    old_Eint_array(k,6) = 0;
    new_Eint_array(k,6) = 0;
end
if neighbourmask(2,2) ~= nospace
    k = neighbourmask(2,2);
    old_Eint_array(k,8) = 0;
    new_Eint_array(k,8) = 0;
end
if (neighbourmask(1,1) ~= nospace) && (neighbourmask(1,2) ~= nospace)
    k = neighbourmask(1,2);  %%correct 2x2 values on the diagonal of Eint_array ...
    old_Eint_array(k,2) = 0;
    new_Eint_array(k,2) = 0;
    k = neighbourmask(1,1);
    old_Eint_array(k,4) = 0;
    new_Eint_array(k,4) = 0;
end
if (neighbourmask(2,1) ~= nospace) && (neighbourmask(2,2) ~= nospace)
    k = neighbourmask(2,2);
    old_Eint_array(k,6) = 0;
    new_Eint_array(k,6) = 0;
    k = neighbourmask(2,1);
    old_Eint_array(k,8) = 0;
    new_Eint_array(k,8) = 0;
end
if (neighbourmask(1,1) ~= nospace) && (neighbourmask(2,1) ~= nospace)
    k = neighbourmask(2,1); %%... and 4x2 values on the offdiagonal
    old_Eint_array(k,2) = 0;
    new_Eint_array(k,2) = 0;
    k = neighbourmask(1,1);
    old_Eint_array(k,6) = 0;
    new_Eint_array(k,6) = 0;
end
if (neighbourmask(1,2) ~= nospace) && (neighbourmask(2,1) ~= nospace)
    k = neighbourmask(2,1);
    old_Eint_array(k,4) = 0;
    new_Eint_array(k,4) = 0;
    k = neighbourmask(1,2);
    old_Eint_array(k,6) = 0;
    new_Eint_array(k,6) = 0;
end
if (neighbourmask(1,1) ~= nospace) && (neighbourmask(2,2) ~= nospace)
    k = neighbourmask(2,2);
    old_Eint_array(k,2) = 0;
    new_Eint_array(k,2) = 0;
    k = neighbourmask(1,1);
    old_Eint_array(k,8) = 0;
    new_Eint_array(k,8) = 0;
end
if (neighbourmask(1,2) ~= nospace) && (neighbourmask(2,2) ~= nospace)
    k = neighbourmask(2,2);
    old_Eint_array(k,4) = 0;
    new_Eint_array(k,4) = 0;
    k = neighbourmask(1,2);
    old_Eint_array(k,8) = 0;
    new_Eint_array(k,8) = 0;
end

old_Eint = sum(sum(old_Eint_array)); 
new_Eint = sum(sum(new_Eint_array));


%% Compute new Eext for the cell

old_Eext_array = zeros(2,2);
new_Eext_array = zeros(2,2);

if neighbourmask(1,1) ~= nospace
    old_Eext_array(1,1) = 0.25 * my * (deltaX*deltaY) * (G_old(1,1)+G_old(1,2)+G_old(2,1)+G_old(2,2));
    new_Eext_array(1,1) = 0.25 * my * (deltaX*deltaY) * (G(1,1)+G(1,2)+G(2,1)+G(2,2));
end
if neighbourmask(1,2) ~= nospace
    old_Eext_array(1,2) = 0.25 * my * (deltaX*deltaY) * (G_old(1,2)+G_old(1,3)+G_old(2,2)+G_old(2,3));
    new_Eext_array(1,2) = 0.25 * my * (deltaX*deltaY) * (G(1,2)+G(1,3)+G(2,2)+G(2,3));
end
if neighbourmask(2,1) ~= nospace
    old_Eext_array(2,1) = 0.25 * my * (deltaX*deltaY) * (G_old(2,1)+G_old(2,2)+G_old(3,1)+G_old(3,2));
    new_Eext_array(2,1) = 0.25 * my * (deltaX*deltaY) * (G(2,1)+G(2,2)+G(3,1)+G(3,2));
end
if neighbourmask(2,2) ~= nospace
    old_Eext_array(2,2) = 0.25 * my * (deltaX*deltaY) * (G_old(2,2)+G_old(2,3)+G_old(3,2)+G_old(3,3));
    new_Eext_array(2,2) = 0.25 * my * (deltaX*deltaY) * (G(2,2)+G(2,3)+G(3,2)+G(3,3));
end

old_Eext = Ha * sum(sum(old_Eext_array));
new_Eext = Ha * sum(sum(new_Eext_array));


%% Compute new Ekin for the cell

old_Ekin_array = 0.5 * my * Lambda * (old_gradGX.*old_gradGX + old_gradGY.*old_gradGY)*(deltaX*deltaY); 
new_Ekin_array = 0.5 * my * Lambda * (new_gradGX.*new_gradGX + new_gradGY.*new_gradGY)*(deltaX*deltaY); 

old_Ekin = sum(sum(old_Ekin_array));
new_Ekin = sum(sum(new_Ekin_array));


%% New energies

dEint = new_Eint - old_Eint;
dEext = new_Eext - old_Eext;
dEkin = new_Ekin - old_Ekin;

dE = dEint + dEext + dEkin;


end






