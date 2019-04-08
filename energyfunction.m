%% This script computes the Energy of the starting G
%%Istead of looping through the chip, this algorithm usus matrix
%%manipulation to compute the energy (which is much fatser than looping for matlab) 

function [E,gradGY_vec,gradGX_vec] = energyfunction(G,coord_array,Nconst,Ha,Lambda,geometry)


my = 4*pi*1e-07;                   %%a physical constant!                

gridpointX = geometry.gridpointX;
gridpointY = geometry.gridpointY;
deltaX = geometry.deltaX;
deltaY = geometry.deltaY;
nospace = geometry.nospace;

coord_vecdummy = coord_array';
coord_vec = coord_vecdummy(:);


%% the 4 arrays G1, G2, G3, and G4 are needed for the following calculations
%  for a certain index (i,j) they contain the 4 gridpoints that belong to a cell

Gdummy = G;
Gdummy(gridpointY,:) = [];
Gdummy(:,gridpointX) = [];
Gdummy(coord_array==nospace) = 0;
G1 = Gdummy;    %%for the gridpoint in the upper left corner of a cell
Gdummy = G;
Gdummy(gridpointY,:) = [];
Gdummy(:,1) = [];
Gdummy(coord_array==nospace) = 0;
G2 = Gdummy;    %%for the gridpoint in the upper right corner of a cell
Gdummy = G;
Gdummy(1,:) = [];
Gdummy(:,gridpointX) = [];
Gdummy(coord_array==nospace) = 0;
G3 = Gdummy;    %%for the gridpoint in the left down corner of a cell
Gdummy = G;
Gdummy(1,:) = [];
Gdummy(:,1) = [];
Gdummy(coord_array==nospace) = 0;
G4 = Gdummy;    %%for the gridpoint in the right down corner of a cell


%% Firstly we compute the gradients of G in the y and x direction by bilinear interpolation

gradGY = (0.5/deltaY) * ((G3-G1)+(G4-G2));
gradGX = (0.5/deltaX) * ((G2-G1)+(G4-G3));

gradGYdummy = gradGY';
gradGXdummy = gradGX';
gradGY_vec = gradGYdummy(:);
gradGX_vec = gradGXdummy(:);
gradGY_vec = gradGY_vec(coord_vec ~= 0);
gradGX_vec = gradGX_vec(coord_vec ~= 0);


%% 1) Compuatuion of Eint: interaction of induced currents and the field created by themselves 
Eint_array = 0.5 * my * ((gradGY_vec*gradGY_vec'+gradGX_vec*gradGX_vec').*Nconst);
Eint = sum(sum(Eint_array));

%% 2) Computation of Eext: interaction of currents with the external field Ha    
Eext_array = 0.25 * my * Ha * (deltaX*deltaY) * (G1 + G2 + G3 + G4);
Eext = sum(sum(Eext_array)); 

%% 3) Computation of Ekin: kinetic energy of charge carriers 
Ekin_array = 0.5 * my * Lambda * (deltaX*deltaY) * (gradGX.*gradGX + gradGY.*gradGY); 
Ekin = sum(sum(Ekin_array));


%% And at least, the sum of them all:

E = Eext + Ekin + Eint;

end

