function t = ant_to_pos(t,iPoles,side)

ns    = length(t)-1;
apole = iPoles(1);
ppole = iPoles(2);

t = t(1:ns);
if side == 1        % anterior pole to posterior
    if apole < ppole
        t = t(apole+1:ppole);
    else
        t = [t(apole+1:ns);t(1:ppole)];
    end
elseif side == 2    % posterior to anterior flipped    
    if apole < ppole
        t = [flipud(t(1:apole-1)); flipud(t(ppole:ns))];
    else
        t = flipud(t(ppole:apole-1));
    end
end


end