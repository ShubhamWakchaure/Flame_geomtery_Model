function y = distance(p,q)
    x_diff = p(1)-q(1);
    y_diff = p(2)-q(2);
    y = sqrt( (x_diff)^2 + (y_diff)^2 );
    
end