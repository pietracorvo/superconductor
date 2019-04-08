%% This script runs some preliminary steps
%
% Here we make some geometrical analysis of the chip and implement the bundary conditions: 
% --> "geometrymask" is the G sized mask which notes the conducting parts and holes
% --> "currentmask" is the G sized mask which notes the current boundary conditions
% --> "holemask" is the G sized mask which enumerates the holes in the chip
% --> "coord" is a three column vector, with the cell index, its row and its column index
%     "coord_array" a (gridpointX-1)*(gridpointY-1) sized array with the cellnumbers in it
% --> "neighbourmask" is a 2*2*(gridpointX-1)*(gridpointY-1) array, that is
%       needed for the delta energy method ...
% input: \
% output: G, geometrymask, coord, coord_array, neighbourmask, geometry(a structure with all relevant 
%         gemoetrical properties and codewords), currentmask, holemask and whole_curent, Ha  


%% Some cleaning work

E_vec = [];
clear coord;


%% Some physical constants, some definitions 

my = 4*pi*1e-07;   %%a physical constant ... 

%%physical parameters
thickness = 0.5*1e-06;                    %%thickness of the superconductor [m]
penetration_depth = 35e-09;               %%london penetration depth [m]
Lambda = penetration_depth^2/thickness;   %%characteristical parameter of superconductor [m]
W = 10^(-5);

%%'codewords' for the arrays:
nospace = 0;   %%this is how we note 'there is no conducter but not a hole' in geometrymask 
space = 1;     %%this is how we note 'there is conductor' in geometrymask
hole = -2;     %%this is how we note 'there is a hole' in geometrymask
nohole = 0;    %%this is how we note 'that is not a hole' in holemask, all other numbers in holemask are holenumbers

nocurrent = Inf; %%this is how we note 'here is no current boundary condition' in currentmask and general_currentmask
general_current_1 = 1;  %%these two numbers symbolize the 2 different kinds of borders in general_currentmask
general_current_2 = 2;  %%          - " - 


%% This are some basic parameters

% starting= questdlg('How do you want to start?', 'Getting started', 'input parameters', 'standard parameters', 'standard parameters');
starting = 'standard parameters';
switch starting
    case 'input parameters'
        disp('Do you want to give me the dimensions and the number of gridpoints?');
        dX = input('Length in the X direction [m]: ');
        dY = input('Length in the Y direction [m]: ');
        gridpointX = input('Number of gridpoints in the X direction: ');
        gridpointY = input('Number of gridpoints in the Y direction: ');
    case 'standard parameters'
        % the wires have a thickness if 10 my 
        dX = 1.96*10^(-5); %%this are m                               
        dY = 1.96*10^(-5); %%   - " -
        gridpointX = 50;                      
        gridpointY = 50;
        disp(['Length in the X direction [m]: ' num2str(dX)]);
        disp(['Length in the Y direction [m]: ' num2str(dY)]);
        disp(['Number of gridpoints in the X direction: ' num2str(gridpointX)]);
        disp(['Number of gridpoints in the Y direction: ' num2str(gridpointY)]);
end

deltaX = dX/(gridpointX-1);                    %%here we compute the dimension of one single cell                                
deltaY = dY/(gridpointY-1);                    %%                       -  "  -

starter = 1e-6; %that is to have nonzero starting energy in the magnetic case
G = starter * ones(gridpointY,gridpointX);    %%now make a G with some background noise 


%% Here we get the special geometry

% make_geometry = questdlg('Which geometry you want to use?', 'Geometry input', 'strip', 'halfpiturn', 'halfpiturn', 'none');
make_geometry = 'strip';
switch make_geometry
    case 'strip'
        %%thats the strip in positive x direction
        geometrymask = space * ones(gridpointY,gridpointX);
        geometrymask(1:ceil(gridpointY*0.25),:) = nospace;
        geometrymask(ceil(gridpointY*0.75):gridpointY,:) = nospace;  
    case 'halfpiturn'
        make_halfpiturn;
    case 'island'
        geometrymask = nospace * ones(gridpointY,gridpointX);
        geometrymask(13:38,13:38) = space;
    case 'strip_hole'
        %%thats the strip in positive x direction with a hole in the middle
        geometrymask = space * ones(gridpointY,gridpointX);
        geometrymask(1:ceil(gridpointY*0.25),:) = nospace;
        geometrymask(ceil(gridpointY*0.75):gridpointY,:) = nospace;  
        geometrymask(ceil(gridpointY*0.4):ceil(gridpointY*0.6),ceil(gridpointX*0.4):ceil(gridpointX*0.6)) = nospace; %whole
end

%figure; imagesc(geometrymask); colorbar,


%% This function 'holefinder' gives a mask called 'holemask' with the holes

[ holemask, holenumber ] = holefinder(geometrymask, gridpointX, gridpointY, nospace, nohole, space); 
geometrymask(holemask~=nohole)= hole;  

%%So there are two ways of 'nothing' in geometrymask:
%%. 1) 'geometrymask == space' means conductor regions without holes ('geometrymask == space' is nonconducting regions without holes)
%%. 2) 'geometrymask ~= nospace' means conductor regions with holes (we need this for all calculations...)


%% Here we fit the G with the geometrymask and holemask

G(geometrymask==nospace) = 0;

%%fit the same value of G for holeregions
if (holenumber ~= 0)
for i = 1:1:holenumber
    holeval = sum(G(holemask == 1)) / size(G(holemask == 1),1);
    G(holemask==i) = holeval;
end
end


%% Here we implement coord and coord_array which connect the k index with space coordinates

k = 0;
coord_array = nospace * ones(gridpointY-1,gridpointX-1);
for i = 1:1:gridpointY-1
for j = 1:1:gridpointX-1
    if (geometrymask(i,j) ~= nospace) && (geometrymask(i+1,j) ~= nospace) && ...
       (geometrymask(i,j+1) ~= nospace) && (geometrymask(i+1,j+1) ~= nospace) 
        k = k+1;
        coord(k,1) = k; %#ok<SAGROW>    %%1. row: k for existing cell
        coord(k,2) = i; %#ok<SAGROW>    %%2. row: y coordnate of the cell k
        coord(k,3) = j; %#ok<SAGROW>    %%3. row: x coordnate of the cell k
        coord_array(i,j) = k;           %%this converts space index to cell index
    end
end
end
real_cellnumber = k;    %%the number of counducting cells (is the size of 'coord') 


%% This makes you the characteristical structure 'geometry', it is just to make function calling more handsome ...

%%and now we fit a structure called 'geometry', with all these geometrical information and the codewords
geometry = struct( 'dX',dX, 'dY',dY, 'gridpointX',gridpointX, 'gridpointY',gridpointY, 'thickness', thickness, ...
            'deltaX', deltaX, 'deltaY',deltaY, 'space',space, 'nospace', nospace, 'nohole',nohole, 'hole',hole, ...
            'nocurrent', nocurrent, 'general_current_1',general_current_1, 'general_current_2',general_current_2, ...
            'real_cellnumber', real_cellnumber, 'holenumber', holenumber );          


%% Analyse geometrymask for the neighbourhood of a gridpoint

%%For every grdipoint (i,j) there are 4 surrounding cells
%%when this cells don't exist, they are noted with 'nospace' in 'neighbourmask',
%%all other (natural) numbers in 'neighbourmask', are the numbers of the cells

neighbourmask = nospace * ones(2,2,gridpointY,gridpointX);

neighbourmask(1,1,2:gridpointY,2:gridpointX) = coord_array;
neighbourmask(1,2,2:gridpointY,1:gridpointX-1) = coord_array;
neighbourmask(2,1,1:gridpointY-1,2:gridpointX) = coord_array;
neighbourmask(2,2,1:gridpointY-1,1:gridpointX-1) = coord_array;


%% Here we get boundary conditions input: Ha and the external current

%boundary_condition_1 = questdlg ('Is there any external current?', 'Boundary conditions: External currents', 'Yes', 'No', 'Standard current', 'Standard current');
boundary_condition_1 = 'Standard current';
switch boundary_condition_1
    case 'Yes'
        whole_current = input('How much current is comming in [A]: ');  
    case 'No'
        whole_current = 0;
        disp(['No external currents applied.' '']);
    case 'Standard current' 
        whole_current = 1;
        disp(['Standard current: ' num2str(whole_current) ' A']);
end

%%this function looks for the location of currents and outputs general_currentmask
[general_currentmask,~] = currents(geometrymask,holemask,geometry);

%%make a currentmask from general_currentmask and writes the value of the currents in it
current_2 = whole_current/2+starter;
current_1 = (-1) * whole_current/2+starter;
currentmask = general_currentmask;
currentmask(currentmask==general_current_2) = current_2;
currentmask(currentmask==general_current_1) = current_1;

%%Write the boundary conditions in currrentmask in the G matrix
G(currentmask==current_2) = current_2;
G(currentmask==current_1) = current_1;


%boundary_condition_2 = questdlg ('Is there any external field?', 'Boundary conditions: External field Ha ', 'Yes', 'No','Standard field', 'No');         
boundary_condition_2 = 'Standard field';
switch boundary_condition_2
    case 'Yes'
        Ha = input('How strong is the external magnetic field [A/m]: ');
        case 'No'
        disp(['No external field applied.' '']);
        Ha = 0;
    case 'Standard field'
        Ha = 10^5;
        disp(['Standard field: ' num2str(Ha) ' A/m']);
end


%% The numerical soultions of N

tic
Nconst = zeros(real_cellnumber,real_cellnumber);
for k = 1:1:real_cellnumber
    for k2 = 1:1:real_cellnumber 
        Nconst(k,k2) = N(k,k2,geometry,coord); 
    end
end              
time.t_Nconst = toc;


%%

disp('ready')
