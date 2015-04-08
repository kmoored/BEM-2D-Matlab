function [ x, y ] = rotatePts( x_0, y_0, theta)
%rotate: Rotates a set of points a specified angle in-plane.
    x = x_0 .* cos(theta) - y_0 .* sin(theta);
    y = x_0 .* sin(theta) + y_0 .* cos(theta);
end

