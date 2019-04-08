%% Change The values in G and the gradient vectors
%%After applied variation we have to make a new G_slice
%%and new gradient vectors, used for the energy calculation by the 'delta' method


%% Make new gradient vectors(with energyfunction style ...)

coord_vecdummy = coord_array';
coord_vec = coord_vecdummy(:);

Gdummy = G;
Gdummy(gridpointY,:) = [];
Gdummy(:,gridpointX) = [];
Gdummy(coord_array==geometry.nospace) = 0;
G1 = Gdummy;    %%for the gridpoint in the upper left corner of a cell
Gdummy = G;
Gdummy(gridpointY,:) = [];
Gdummy(:,1) = [];
Gdummy(coord_array==geometry.nospace) = 0;
G2 = Gdummy;    %%for the gridpoint in the upper right corner of a cell
Gdummy = G;
Gdummy(1,:) = [];
Gdummy(:,gridpointX) = [];
Gdummy(coord_array==geometry.nospace) = 0;
G3 = Gdummy;    %%for the gridpoint in the left down corner of a cell
Gdummy = G;
Gdummy(1,:) = [];
Gdummy(:,1) = [];
Gdummy(coord_array==geometry.nospace) = 0;
G4 = Gdummy;    %%for the gridpoint in the right down corner of a cell

gradGY = (0.5/geometry.deltaY) * ((G3-G1)+(G4-G2));
gradGX = (0.5/geometry.deltaX) * ((G2-G1)+(G4-G3));

gradGYdummy = gradGY';
gradGXdummy = gradGX';
gradGY_vec = gradGYdummy(:);
gradGX_vec = gradGXdummy(:);
gradGY_vec = gradGY_vec(coord_vec ~= 0);
gradGX_vec = gradGX_vec(coord_vec ~= 0);


%% Slice the G

%%for an index (i,j) in G_slice there are small 3x3 block matrices
%%which represent the neighbourhood of this gridpoint

G_slice = zeros(3,3,gridpointY,gridpointX);
for i2 = 1:1:gridpointY
for j2 = 1:1:gridpointX

    G_slice(2,2,i2,j2) = G(i2,j2); %%in the middle of the 3x3 array we plug the same number as G
    
    %%and now we look if the neighbourhood of this gridpoint exists and
    %%copy their values in the small 3x3 blocks in G_slice
    if neighbourmask(1,1,i2,j2) ~= geometry.nospace
        G_slice(1,1,i2,j2) = G(i2-1,j2-1);
        G_slice(2,1,i2,j2) = G(i2,j2-1);
        G_slice(1,2,i2,j2) = G(i2-1,j2);
    end
    if neighbourmask(1,2,i2,j2) ~= geometry.nospace
        G_slice(1,3,i2,j2) = G(i2-1,j2+1);
        G_slice(1,2,i2,j2) = G(i2-1,j2);
        G_slice(2,3,i2,j2) = G(i2,j2+1);
    end
    if neighbourmask(2,1,i2,j2) ~= geometry.nospace
        G_slice(3,1,i2,j2) = G(i2+1,j2-1);
        G_slice(2,1,i2,j2) = G(i2,j2-1);
        G_slice(3,2,i2,j2) = G(i2+1,j2);
    end
    if neighbourmask(2,2,i2,j2) ~= geometry.nospace
        G_slice(3,3,i2,j2) = G(i2+1,j2+1);
        G_slice(2,3,i2,j2) = G(i2,j2+1);
        G_slice(3,2,i2,j2) = G(i2+1,j2);
    end
    
end
end


