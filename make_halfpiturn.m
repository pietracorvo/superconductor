%% Sharp pi turn geometry

%%input: gridpointX, gridpointY
%%output: geometrymask

nospace = 0;
geometrymask = ones(gridpointY,gridpointX);

whole_start_X = 1; 
whole_start_Y = ceil(gridpointY/2);
whole_end_X = floor(gridpointX/2);
whole_end_Y = gridpointY;

for i = whole_start_Y+3:1:whole_end_Y
    for j = whole_start_X:1:whole_end_X-2
        geometrymask(i,j) = nospace;
    end
end

geometrymask(1,:) = nospace;
geometrymask(:,gridpointX) = nospace;