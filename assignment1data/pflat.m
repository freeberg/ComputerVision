function [y]=pflat(x);

    for kol=1:size(x,2)
        y(:,kol) = x(:,kol)./x(end,kol);
    end

    
