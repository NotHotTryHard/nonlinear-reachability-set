function [t_end, X_tmp, Y_tmp] = find_cross_x1(alpha, x0, T1, sg) %ищет t1 и t2

    eps = 10^-5;

    hold on
    grid on

    t0 = 0;
    t_end = 0;
    X_tmp = [];
    Y_tmp = [];
    while 1

        if sg == -1
            options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @event_cross_x1);
            [t, x] = ode45(@(t, x) S_minus(t, x, alpha), linspace(t0, T1, T1 * 100), x0, options);
        elseif sg == 1
            options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @event_cross_x1);
            [t, x] = ode45(@(t, x) S_plus(t, x, alpha), linspace(t0, T1, T1 * 100), x0, options);
        end
        X_tmp = [X_tmp; x(:, 1)];
        Y_tmp = [Y_tmp; x(:, 2)];


        if abs(t(end) - T1) < eps
            break;
        end
        if t_end == 0
            t_end = t(end);
        end
        % plot(x(end, 1), 0, '*')
        x0 = x(end, :);
        x0(2) = 0;
        x0(4) = 0;
        t0 = t(end);
        sg = -sg;
        break
    end
end

function [value, isterminal, direction] = event_cross_x1(t, x)
    value = [x(2)];
    isterminal = [1];
    direction = [0];
end
