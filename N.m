%% Numerical solution of the double surface integral N(k,k2)
%
%You shouldn't be here, that's just some kind of ominous mathematical things...


function [valueN] = N(k,kp,geometry,coord) 

deltaX = geometry.deltaX;
deltaY = geometry.deltaY;
dX = geometry.dX;
dY = geometry.dY;

%%that makes 2 spatial indices in x and y for a cellindex k
    y = coord(k,2);
    x = coord(k,3);
    yp = coord(kp,2);
    xp = coord(kp,3);
    
%%this here converts the 'space indices' y,x,yp,xp to 'space coordinates'    
    y_vec = deltaY/2:deltaY:dY-deltaY/2;
    x_vec = deltaX/2:deltaX:dX-deltaX/2;
    y = y_vec(y);
    x = x_vec(x);
    yp = y_vec(yp);
    xp = x_vec(xp);
    dy = deltaY;
    dx = deltaX;
    
    Nkkp(1)  =  T(xp+dx,x+dx,yp+dy,y+dy); %%you see what I said: you don't want to take a look at that ...
    Nkkp(2)  = -T(xp+dx,x+dx,yp+dy,y);
    Nkkp(3)  = -T(xp+dx,x,yp+dy,y+dy);
    Nkkp(4)  =  T(xp+dx,x,yp+dy,y);
    Nkkp(5)  = -T(xp,x+dx,yp+dy,y+dy);
    Nkkp(6)  =  T(xp,x+dx,yp+dy,y);
    Nkkp(7)  =  T(xp,x,yp+dy,y+dy);
    Nkkp(8)  = -T(xp,x,yp+dy,y);
    Nkkp(9)  = -T(xp+dx,x+dx,yp,y+dy);
    Nkkp(10) =  T(xp+dx,x+dx,yp,y);
    Nkkp(11) =  T(xp+dx,x,yp,y+dy);
    Nkkp(12) = -T(xp+dx,x,yp,y);
    Nkkp(13) =  T(xp,x+dx,yp,y+dy);
    Nkkp(14) = -T(xp,x+dx,yp,y);
    Nkkp(15) = -T(xp,x,yp,y+dy);
    Nkkp(16) =  T(xp,x,yp,y);

    valueN = sum(Nkkp);
end

function valueT = T(xp,x,yp,y)
    % first term
    ft = -(1/6)*((x-xp)^2+(y-yp)^2)^(3/2);
    
    % prefactor
    pf = 0.5*(x-xp)*(y-yp);
    
    % first log
    flog = (y-yp)*log(x-xp+sqrt((x-xp)^2+(y-yp)^2));
        
    % second log
    slog = (x-xp)*log(y-yp+sqrt((x-xp)^2+(y-yp)^2));
    
    if isnan(flog) || isinf(flog)
        flog = 0;
    end
    if isnan(slog) || isinf(slog)
        slog = 0;
    end
    
    valueT = ft+pf*(flog+slog);
end
