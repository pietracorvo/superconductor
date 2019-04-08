
%% This programm looks for the borders of your chip and identificates the current boundary condition
%
% This algorithm will just be able to find the borders of your geometry, if there is at least one 
% row or column between it and the chip border parallel to the conductor.(Thatmeans: all conductor, 
% which starts at the border, is for in or outcomming current)
% Currents are counted positive from left to right or from up to down, (like the coordinates) 
% The resolution must be high enough, a straight borderline section must be at leat made of three gridpoints 
% (otherwise corners of the borderlines won't be detected!).  


function [general_currentmask,length_currentzone] = currents(geometrymask,holemask,geometry)

border = -1;

gridpointX = geometry.gridpointX;
gridpointY = geometry.gridpointY;
deltaX = geometry.deltaX;
deltaY = geometry.deltaY;
nospace = geometry.nospace;
space = geometry.space;
nohole = geometry.nohole;
nocurrent = geometry.nocurrent;
general_current_1 = geometry.general_current_1;
general_current_2 = geometry.general_current_2;


%% Initiate an array

general_currentmask = ones(gridpointY,gridpointX) * nocurrent;
        

%% First detect borders and note them in general_currentmask 

for i = 2:1:gridpointY-1
for j = 1:1:gridpointX
%%for the row
    if (geometrymask(i-1,j) == nospace) && (geometrymask(i,j) == space) 
        general_currentmask(i,j) = border;     
    end
    if (geometrymask(i+1,j) == nospace) && (geometrymask(i,j) == space)  
        general_currentmask(i,j) = border;     
    end
end
end

for i = 1:1:gridpointY
for j = 2:1:gridpointX-1
%%for the column
    if (geometrymask(i,j-1) == nospace) && (geometrymask(i,j) == space) 
        general_currentmask(i,j) = border;     
    end
    if (geometrymask(i,j+1) == nospace) && (geometrymask(i,j) == space) 
        general_currentmask(i,j) = border;     
    end
end
end

%%sometimes 'corners' in the boundary condition are not seen by the above iteration, we fix that here ...
%%(but this algorithm won't work with curved borders! minimum length of straight sections: 3 gridpoints)
cornercount = 0;
for i = 3:1:gridpointY-2
for j = 3:1:gridpointX-2
    %%there are 4 possible 'corner cases'
    if general_currentmask(i+1,j)==general_currentmask(i,j+1) ...
        && general_currentmask(i,j)==nocurrent ...
        && general_currentmask(i+1,j)==border ...
        && general_currentmask(i,j+1)==border ...
        && general_currentmask(i+2,j)==border ...
        && general_currentmask(i,j+2)==border
            general_currentmask(i,j) = general_currentmask(i,j+1);
            cornercount = cornercount + 1;
    end
    if general_currentmask(i-1,j)==general_currentmask(i,j-1) ...
        && general_currentmask(i,j)==nocurrent ...
        && general_currentmask(i-1,j)==border ...
        && general_currentmask(i,j-1)==border ...
        && general_currentmask(i-2,j)==border ...
        && general_currentmask(i,j-2)==border
            general_currentmask(i,j) = general_currentmask(i,j-1);
            cornercount = cornercount + 1;
    end
    if general_currentmask(i-1,j)==general_currentmask(i,j+1) ...
        && general_currentmask(i,j)==nocurrent ...    
        && general_currentmask(i-1,j)==border ...
        && general_currentmask(i,j+1)==border ...
        && general_currentmask(i-2,j)==border ...
        && general_currentmask(i,j+2)==border
            general_currentmask(i,j) = general_currentmask(i,j+1);
            cornercount = cornercount + 1;
    end
    if general_currentmask(i,j-1)==general_currentmask(i+1,j) ...
        && general_currentmask(i,j)==nocurrent ...    
        && general_currentmask(i,j-1)==border ...
        && general_currentmask(i+1,j)==border ...
        && general_currentmask(i,j-2)==border ...
        && general_currentmask(i+2,j)==border
            general_currentmask(i,j) = general_currentmask(i+1,j);
            cornercount = cornercount + 1;
    end
end
end


%% Now we note the two different borders

%For the left and right border
for index = [1 gridpointX]
for i = 2:1:gridpointY-1
    if general_currentmask(i,index) == border
        %%identificate class of border
        if (geometrymask(i-1,index) == nospace) && (geometrymask(i,index) == space) && (general_currentmask(i,index) == border) 
            general_current = general_current_1;
        end
        if (geometrymask(i+1,index) == nospace) && (geometrymask(i,index) == space) && (general_currentmask(i,index) == border)   
            general_current = general_current_2;
        end
        %%follow the borderline and apply values
        i2 = i;
        if index == 1; j2 = 2; elseif index == gridpointX; j2 = gridpointX-1; end
        general_currentmask(i,index) = general_current;
        general_currentmask(i,j2) = general_current;
        while (i2 ~= 1) && (j2 ~= 1) && (i2 ~= gridpointY) && (j2~= gridpointX)
            if general_currentmask(i2-1,j2) == border
                general_currentmask(i2-1,j2) = general_current;
                i2 = i2-1;
            end
            if general_currentmask(i2+1,j2) == border
                general_currentmask(i2+1,j2) = general_current;
                i2 = i2+1;
            end
            if general_currentmask(i2,j2-1) == border
                general_currentmask(i2,j2-1) = general_current;
                j2 = j2-1;
            end
            if general_currentmask(i2,j2+1) == border
                general_currentmask(i2,j2+1) = general_current;
                j2 = j2+1;
            end
        end
        
    end
end
end
%%For the upper and lower border
for index = [1 gridpointY]
for j = 2:1:gridpointX-1
    if general_currentmask(index,j) == border
        %%identificate class of border
        if (geometrymask(index,j-1) == nospace) && (geometrymask(index,j) ~= nospace) && (general_currentmask(index,j) == border)  
            general_current = general_current_2; 
        end
        if (geometrymask(index,j+1) == nospace) && (geometrymask(index,j) ~= nospace) && (general_currentmask(index,j) == border)  
            general_current = general_current_1;
        end
        %%follow the borderline and apply values
        j2 = j;
        if index == 1; i2 = 2; elseif index == gridpointY; i2 = gridpointY-1; end
        general_currentmask(index,j) = general_current;
        general_currentmask(i2,j) = general_current;
        while (i2 ~= 1) && (j2 ~= 1) && (i2 ~= gridpointY) && (j2~= gridpointX)
            if general_currentmask(i2-1,j2) == border
                general_currentmask(i2-1,j2) = general_current;
                i2 = i2-1;
            end
            if general_currentmask(i2+1,j2) == border
                general_currentmask(i2+1,j2) = general_current;
                i2 = i2+1;
            end
            if general_currentmask(i2,j2-1) == border
                general_currentmask(i2,j2-1) = general_current;
                j2 = j2-1;
            end
            if general_currentmask(i2,j2+1) == border
                general_currentmask(i2,j2+1) = general_current;
                j2 = j2+1;
            end
        end
        
    end
end
end

%% Here we delete the borders of the holes, they are no current boundary conditions ... 

index = [-1 -1; -1 0; -1 1;0 -1; 0 0; 0 1; 1 -1; 1 0; 1 1];
for i = 2:1:gridpointY-1
for j = 2:1:gridpointX-1
    for k = 1:1:9
        if (general_currentmask(i,j)==-1) && (holemask(i+index(k,1),j+index(k,2))~=nohole)
            general_currentmask(i,j) = nocurrent; 
        end
    end
end
end


%% And check the length of the incomming currentzone

length_currentzone = 0;
flag = 0;
for i = 1:1:gridpointY   %%For the left border ...
    if general_currentmask(i,1)==general_current_1
        flag = 1;
        begin = i;
    end
    if (general_currentmask(i,1)==general_current_2) && (flag==1)
        flag = 0;
        length_currentzone = length_currentzone + deltaY*(i-begin);
    end
end
flag = 0;
for j = 1:1:gridpointX   %%... and the upper border.
    if general_currentmask(1,j)==general_current_2
        flag = 1;
        begin = j;
    end
    if (general_currentmask(1,j)==general_current_1) && (flag==1)
        flag = 0;
        length_currentzone = length_currentzone + deltaX*(j-begin);
    end    
end


%% Island geometries (but not holes!) have borders which are not detected from the above algorithm 
%%(because they are not connected to the borders of the chip), we give 'closed' 
%%borderlines (but not hole borders) the boundary condition general_current_1

general_currentmask(general_currentmask == border) = general_current_1; % hier viellciht eine hüpflösung wie oben!!!


%% Check if all -1's are gone

if any(general_currentmask==-1)
   disp('Something with borderdetection went wrong!'); 
end


end

