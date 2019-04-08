%% Here starts one iteration
%
%

%% These parameters must be inititalized before we go into the while loop:

disp('Starting iteration:')

gridpointX = geometry.gridpointX;
gridpointY = geometry.gridpointY;

%%Calculate first energy
t_energy = tic;
[E,~,~] = energyfunction(G,coord_array,Nconst,Ha,Lambda,geometry);   %%calculate the whole energy the first time, later we will just modify the associated tensors ...
time.t_energyfunction = toc(t_energy);
dG_start = 0.01;
dG = dG_start;

%%Initialize some parameters
rescalingflag = 1;
Eneu = E;                                 
E_varP = E * ones(gridpointY,gridpointX);   %%in this arrays write the effect of variation at that point
E_varN = E * ones(gridpointY,gridpointX);   %%                         - " -
E_var_holeP = E * ones(geometry.holenumber);    %%in this arrays we write the effect of variation of holes
E_var_holeN = E * ones(geometry.holenumber);    %%                          - " -


%% And that is the loop

t_dummy = tic;
% Major loop for the iteration steps with stop-conditions
while (iteration<itpar.maxiteration) && ( (abs(E-Eneu)>E*itpar.epsilonfactor_E) || (rescalingflag == 1) )  
    
    iteration = iteration + 1;      
    rescalingflag = 0;              %%this flag supresses the breakdown of iteration after rescaling dG 
    slice_G;                        %%'slice' G and make new gradient vectors (called gradGY_vec and gradGX_vec)
    E = Eneu;                       %%the new energy from the last step, is the energy to compare with for the following  iteration
    
    E_vec(iteration) = E; %#ok<SAGROW> %%this vector is just to watch the convergence behavior of the energy
    
    t_loopstep = tic;
    %%These two loops test the variation of every node of the chip 
    for j = 1:1:gridpointX  
    for i = 1:1:gridpointY    
        if  (currentmask(i,j) == geometry.nocurrent) && (geometrymask(i,j) == geometry.space)  %%avoids variation on non conductiing regions, current boundary cond. and holeregions                 
            
            %Compute the new energies, with sliced arrays ...
            %(G_slice and Nconst_dummy are the parts of the matrices G and Nconst which are needed for the energy calculation)
            G_dummy = G_slice(:,:,i,j);
            Nconst_dummy = [ neighbourmask(1,1,i,j) neighbourmask(1,2,i,j) neighbourmask(2,1,i,j) neighbourmask(2,2,i,j) ];
            Nconst_dummy(Nconst_dummy<1) = 1; 
            
            dEneuP = delta_energyfunction(G_dummy,dG,Ha,Lambda,Nconst(:,Nconst_dummy),gradGY_vec,gradGX_vec,neighbourmask(:,:,i,j),geometry);  %%... for + dG
            dEneuN = delta_energyfunction(G_dummy,-dG,Ha,Lambda,Nconst(:,Nconst_dummy),gradGY_vec,gradGX_vec,neighbourmask(:,:,i,j),geometry); %%... and - dG 
            
            %... and write the E + dE for the variation of every node in arrays. 
            E_varP(i,j) =  E + dEneuP;
            E_varN(i,j) =  E + dEneuN;
        
        end
    end
    end
    
    %%We also want to calculate the energy for the variation of holeregions, we do this with 'energyfunction',
    %%that's a little bit easier (BUT: be aware of possible rounding mistakes, compared to the delta method!)
    if (geometry.holenumber~=0)
        for i = 1:1:geometry.holenumber
            G_hole_P = G;
            G_hole_N = G;
            G_hole_P(holemask==i) = G(holemask==i) + dG;
            G_hole_N(holemask==i) = G(holemask==i) - dG;
            [E_var_holeP(i),~,~] = energyfunction(G_hole_P,coord_array,Nconst,Ha,Lambda,geometry);
            [E_var_holeN(i),~,~] = energyfunction(G_hole_N,coord_array,Nconst,Ha,Lambda,geometry);
        end
        EminP_hole = min(E_var_holeP(:));
        EminN_hole = min(E_var_holeN(:));
    else
        %%In case of no holes in the chip, just give them a big random value (like 1 ...)
        EminP_hole = 1;
        EminN_hole = 1;
    end
    time.t_loopstep = toc(t_loopstep);
    
    %%Now we look for the smallest energy and apply the variation at this point
    EminP = min(E_varP(:));
    EminN = min(E_varN(:)); 
    if  (EminN < E) || (EminP < E) || (EminP_hole < E) || (EminN_hole < E)  
        %%that means we have minimized the energy with a variation or a holevariation
        if (EminP < EminN) && (EminP <= EminP_hole) && (EminP <= EminN_hole) 
            %%EminP is minimizing the energy
            [miny,minx] = find(E_varP==EminP); 
            miny = miny(1); minx = minx(1);
            factor = 1;
            Eneu = EminP;
            if iteration~=1 G(miny(1),minx(1)) = G(miny(1),minx(1)) + dG; end %#ok<SEPEX,SAGROW>
        elseif (EminN <= EminP) && (EminN <= EminP_hole) && (EminN <= EminN_hole) 
            %%EminN is minimizing the energy
            [miny,minx] = find(E_varN==EminN);
             miny = miny(1); minx = minx(1);
            factor = -1;
            Eneu = EminN;
            G(miny(1),minx(1)) = G(miny(1),minx(1)) - dG; %#ok<SAGROW>
        elseif (EminP_hole < EminN_hole) && (EminP_hole < EminP) && (EminP_hole < EminN)
            %%EminP_hole is minimizing the energy
            [holeindex] = find(E_var_holeP==EminP_hole);
            miny = holeindex*(-1); minx = holeindex*(-1);
            factor = 1;
            Eneu = EminP_hole;
            G(holemask==holeindex(1)) = G(holemask==holeindex(1)) + dG; %#ok<SAGROW>
        elseif (EminN_hole <= EminP_hole) && (EminN_hole < EminP) && (EminN_hole < EminN)
            %%EminN_hole is minimizing the energy
            [holeindex] = find(E_var_holeN==EminN_hole);
            miny = holeindex*(-1); minx = holeindex*(-1);
            factor = -1;
            Eneu = EminN_hole;
            G(holemask==holeindex(1)) = G(holemask==holeindex(1)) - dG; %#ok<SAGROW>        
        else disp('Oh no, something went wrong. Ask the admin ...');
        end
    else
        %%that means no variation minimizes the energy, so we have to rescale dG (if it ins't already too small)
        Eneu = E; 
        miny = 0; minx = 0; factor = 0;
        if (dG > dG_start * itpar.epsilonfactor_dG)
            dG = dG * 0.5;      %%if no variation minimizes the energy, we choose a smaller dG
            rescalingflag = 1;  %%that flag prevents termination of iteration due to the epsilon criterium 
            disp('rescale dG ...');
        end
    end
    %%and for the iteration history (holevariations are noted with holenumber*(-1) ... )
    history(iteration,1) = miny; history(iteration,2) = minx; history(iteration,3) = factor*dG;   %#ok<SAGROW>
    
    %%That is just for some output 
    disp(['Step: ' num2str(iteration) '   Energy: ' num2str(Eneu) '      [varaited node: ('...
          num2str(miny(1)) ',' num2str(minx(1)) ') with: ' num2str(factor) '*' num2str(dG) ']']);

    time.t_whole_iteration = time.t_whole_iteration + toc(t_dummy);
      
end



%% So there are 3 cases of 'convergece'

if iteration == itpar.maxiteration
    disp(['Iteration steps (' num2str(iteration) ') exhausted!   [Energy: ' num2str(Eneu) ']']); 
elseif dG <= dG_start * itpar.epsilonfactor_dG
    disp(['Convergence in dG after ' num2str(iteration) ' steps   [Energy: ' num2str(Eneu) ']']);
else
    disp(['Convergence after ' num2str(iteration) ' steps   [Energy: ' num2str(Eneu) ']']);
end


