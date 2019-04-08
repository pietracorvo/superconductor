%% MovieMaker.

[Y, X] = ndgrid(1:size(G,1), 1:size(G,2));
fId = figure;
fId.Position = [0 0 1100 800];

G_old = zeros(size(G,1),size(G,2)); 

%history(:,3) = 10^5*history(:,3); % for fields


%% Set up the movie.
writerObj = VideoWriter('pihalf_current.avi'); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj); 

phi = 0;
picnumb = 50;

for i=1:50000 %size(history,1)
    
    if history(i,1)~=0
   
        G_old(history(i,1),history(i,2)) = G_old(history(i,1),history(i,2)) + history(i,3);
    
        if mod(i,picnumb)==0
            %picnumb = 1.01*picnumb;
        
            %mask = zeros(size(G_old,1),size(G_old,2));
            %mask(history(i,1),history(i,2)) = 1;
            %mask = mask(:);
            %color = repmat([0 0 1],size(mask,1),1);
            %color(mask==1,:) = [1 0 0];
            %radius = repmat(10,size(mask,1),1);
            %radius(mask==1,:) = 20;

            figure(fId); % Makes sure you use your desired frame.
            scatter3(X(:), Y(:), G_old(:)); % radius, color,'filled'
            title(['iteration: ' num2str(i) '     \DeltaG: ' num2str(abs(history(i,3)))]);
            
            phi=phi+1;
            view(phi,20);

            frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
            writeVideo(writerObj, frame);
        
        end
    end
 
end

hold off
close(writerObj); % Saves the movie.
