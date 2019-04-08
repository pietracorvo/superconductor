%% This script detects holes in geometrymask and notes them with a number in 'holemask'

function [holemask, holenumber] = holefinder(geometrymask, gridpointX, gridpointY, nospace, nohole, space)


holemask = (-1) * ones(gridpointY,gridpointX); 
%%. -1 is the temporary codeword for 'hey, there could be a hole!'
%%After the next operations all -1's will vanish, and at the end
%%there should just be 'nohole' and numbers from 1 to the number of holes,
%%which mark every hole. 


%% First: no-conductor regions, which are connected to the border, are no holes ... 
%%We identificate them like this: we go down one border and note regions without
%%geometry, then we go rowwise or columnwise into the chip and check if the
%%gridpoints in the new row/column are connected with region connected to
%%the nohole region detected during the last step.
%%We need to check this for every border, that measn 4 times.

%%left border
for i = 1:1:gridpointY
    if geometrymask(i,1) == nospace
        holemask(i,1) = nohole; 
    end
end
for j = 2:1:gridpointX 
    flag = 1;
    while (flag == 1)
        flag = 0;
        for i = 2:1:gridpointY-1    
            if (holemask(i,j-1) == nohole) && (geometrymask(i,j) == nospace) && (holemask(i,j) == -1)
                holemask(i,j) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i-1,j) == nospace) && (holemask(i-1,j) == -1)
                holemask(i-1,j) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i+1,j) == nospace) && (holemask(i+1,j) == -1)
                holemask(i+1,j) = nohole;
                flag = 1;
            end
        end
    end
end

%%right border
for i = 1:1:gridpointX
    if geometrymask(i,gridpointX) == nospace
        holemask(i,gridpointX) = nohole; 
    end
end
for j = gridpointX-1:-1:1 
    flag = 1;
    while (flag == 1)
        flag = 0;
        for i = 2:1:gridpointY-1    
            if (holemask(i,j+1) == nohole) && (geometrymask(i,j) == nospace) && (holemask(i,j) == -1)
                holemask(i,j) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i-1,j) == nospace) && (holemask(i-1,j) == -1)
                holemask(i-1,j) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i+1,j) == nospace) && (holemask(i+1,j) == -1)
                holemask(i+1,j) = nohole;
                flag = 1;
            end
        end
    end
end

%%upper border
for j = 1:1:gridpointX
    if geometrymask(1,j) == nospace
        holemask(1,j) = nohole; 
    end
end
for i = 2:1:gridpointY 
    flag = 1;
    while (flag == 1)
        flag = 0;
        for j = 2:1:gridpointX-1    
            if (holemask(i-1,j) == nohole) && (geometrymask(i,j) == nospace) && (holemask(i,j) == -1)
                holemask(i,j) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i,j-1) == nospace) && (holemask(i,j-1) == -1)
                holemask(i,j-1) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i,j+1) == nospace) && (holemask(i,j+1) == -1)
                holemask(i,j+1) = nohole;
                flag = 1;
            end
        end
    end
end

%%lower border
for j = 1:1:gridpointX
    if geometrymask(gridpointY,j) == nospace
        holemask(gridpointY,j) = nohole; 
    end
end
for i = gridpointY-1:-1:1 
    flag = 1;
    while (flag == 1)
        flag = 0;
        for j = 2:1:gridpointX-1    
            if (holemask(i+1,j) == nohole) && (geometrymask(i,j) == nospace) && (holemask(i,j) == -1)
                holemask(i,j) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i,j-1) == nospace) && (holemask(i,j-1) == -1)
                holemask(i,j-1) = nohole;
                flag = 1;
            end
            if (holemask(i,j) == nohole) && (geometrymask(i,j+1) == nospace) && (holemask(i,j+1) == -1)
                holemask(i,j+1) = nohole;
                flag = 1;
            end
        end
    end
end


%% Second: conducting regions are also no holes ...

holemask(geometrymask == space) = nohole;


%% Third: Classify compact holeregions 
%%in 'holemask' there are holes where is written '-1' 
%%Now we have to put the connected regions togethrer and note every hole with 
%%a number from 1 to the number of holes.
%%We go trough the chip. If we hit a '-1' we do the same proecdure like
%%before, but now to identificate connected regions of -1. The -1 will be
%%replaced with the number of the hole.

holecount = 0; %%counter for the holes ...

for i = 2:1:gridpointY-1
for j = 2:1:gridpointX-1
if holemask(i,j) == -1;
%%So: if you hit a hole, put all adjacent -1's together and name them with the number 'holecount'
    holecount = holecount + 1;

    holemask (i,j) = holecount;    
    i2 = i;
    j2 = j;
    
    %%left to right:
    flagx = 1;
    while(flagx == 1)
        flagx = 0;
        flagy = 1;
        while(flagy == 1)           
            if (holemask(i2,j2)==holecount) && (holemask(i2-1,j2)==-1) 
                holemask(i2-1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2+1,j2)==-1) 
                holemask(i2+1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2-1)==-1)
                holemask(i2,j2-1) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2+1)==-1)
                holemask(i2,j2+1) = holecount;
                flagx = 1;
                tmpi2 = i2;
            end
            if holemask(i2+1,j2)==nohole
                flagy = 0;
            end
            i2 = i2 + 1;
        end
        flagy = 1;
        while(flagy == 1)           
            if (holemask(i2,j2)==holecount) && (holemask(i2-1,j2)==-1) 
                holemask(i2-1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2+1,j2)==-1) 
                holemask(i2+1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2-1)==-1)
                holemask(i2,j2-1) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2+1)==-1)
                holemask(i2,j2+1) = holecount;
                flagx = 1;
                tmpi2 = i2;
            end
            if holemask(i2-1,j2)==nohole
                flagy = 0;
            end
            i2 = i2 - 1;
        end
        j2 = j2 + 1;
        i2 = tmpi2;
    end

    %%right to left:
    flagx = 1;
    while(flagx == 1)
        flagx = 0;
        flagy = 1;
        while(flagy == 1)
            flagy = 0;
            if (holemask(i2,j2)==holecount) && (holemask(i2-1,j2)==-1) 
                holemask(i2-1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2+1,j2)==-1) 
                holemask(i2+1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2-1)==-1)
                holemask(i2,j2-1) = holecount;
                flagx = 1;
                tmpi2 = i2;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2+1)==-1)
                holemask(i2,j2+1) = holecount;
            end
            if holemask(i2+1,j2)==nohole
                flagy = 0;
            end
            i2 = i2 + 1;
        end
        flagy = 1;
        while(flagy == 1)
            flagy = 0;
            if (holemask(i2,j2)==holecount) && (holemask(i2-1,j2)==-1) 
                holemask(i2-1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2+1,j2)==-1) 
                holemask(i2+1,j2) = holecount;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2-1)==-1)
                holemask(i2,j2-1) = holecount;
                flagx = 1;
                tmpi2 = i2;
            end
            if (holemask(i2,j2)==holecount) && (holemask(i2,j2+1)==-1)
                holemask(i2,j2+1) = holecount;
            end
            if holemask(i2-1,j2)==nohole
                flagy = 0;
            end
            i2 = i2 - 1;
        end
        j2 = j2 - 1;
        i2 = tmpi2;
    end
    
end        
end
end

%% Define number of holes

holenumber = holecount;


end




