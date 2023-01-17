function [x_out,y_out] = rotatePoints(theta, x, y)

    x_out = x*cos(theta) - y*sin(theta);
    y_out = y*cos(theta) + x*sin(theta);

end