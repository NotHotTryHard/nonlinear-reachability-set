function f = S_plus(t, y, alpha)
    f = [y(2);
            -3*y(1)^3 - 2*sin(y(1)) - y(2)^2 + alpha;
            -(y(4)*(-9*y(1)^2 - 2*cos(y(1))));
            -(y(3) - 2*y(2)*y(4))];
end