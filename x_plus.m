function f = x_plus(t, y, alpha)
    f = [y(2);
            -3*y(1)^3 - 2*sin(y(1)) - y(2)^2 + alpha];
end