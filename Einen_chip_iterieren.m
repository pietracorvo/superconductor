%% Main programm: this torso makes the iteration and variation of G
% there are 3 sections: 1) get the input (getinp.m): 
%                             this gives the starting G for the algorithm
%                       2) the main algorithm (iterate_G.m): 
%                             here we variate the values of G and save the variations 
%                             which minimize the energy
%                       3) plotter (plotter.m): 
%                             plots the resulting currents, magnetic field, ...
% 
% These 3 parts are build in a modular way: that means iterate_G.m (which may run for hours and days)
% can always be stopped and restarted, to calculate on if all variables needed (G, geomtry, 
% geometrymask, currentmask, ...) are present in the workspace. Also the script plotter.m can 
% always be executed, if all variables needed are present in the workspace.
% (just stop the calculations after some thousand steps with Ctrl+C, plot with 'plotter' to see if all
% runs fine and then calculate on by starting 'iterate_G')
% (var_clear.m is a script that deletes garbage not needed for the computation wich collets in the
% workspace when killing the iteration ...)


%% With this function we get the input

getinp;   %%this starts the input script, which analyzes the geometry, the boundary conditions, ...


%% Here begins the algorithm

%%Some iteration parameters: determine convergence behaviour
itpar.maxiteration = 100000;  % structure 'itpar.' -- iteration paramters
itpar.epsilonfactor_E = 1e-20;
itpar.epsilonfactor_dG = 1e-30;

iteration = 0;
time.t_whole_iteration = 0;

iterate_G; %%this script does the variation of G


%% Here begins the output

plotter;



