function y = normalize021(x)
    if length(x) == 1 || all(x==0,'all')
        y = x;
        return
    end
    minx = min(x,[],'all');
    maxx = max(x,[],'all');
    
    y    = (x - minx)/(maxx - minx);
end

