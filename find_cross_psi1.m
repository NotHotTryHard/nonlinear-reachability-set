function [Xfin, data, curve1, curve2] = find_cross_psi1(alpha, f, t, T1, sg_start, N, curve1, curve2) %ищет t1 и t2

    eps = 10^-7;
    tspan = linspace(0, t, N+1);
    [tbeg, ybeg] = ode45(@(t, x) f(t, x, alpha), tspan(1:(end-1)), [0,0]');
    Xfin = [];

    data = [];
    for i = 2:N
        if sg_start == -1
            y0 = [ybeg(i, 1); ybeg(i, 2); 1; 0];
        end
        if sg_start == 1
            y0 = [ybeg(i, 1); ybeg(i, 2); -1; 0];
        end

        % plot(x0(1), x0(2), '.b', MarkerSize=10); %точки выхода из кривой - синие 
        tbeg_tmp = tbeg(i);
        sg = sg_start;
        lap = 1;
        while 1
            
            if sg == -1
                options = odeset('RelTol' , 1e-4, 'AbsTol' , 1e-7, 'Events', @events_cross_psi1);
                [t, y] = ode45(@(t, y) S_minus(t, y, alpha), linspace(tbeg_tmp, T1, 100*T1), y0, options); 
            else
                options = odeset('RelTol' , 1e-4, 'AbsTol' , 1e-7, 'Events', @events_cross_psi1);
                [t, y] = ode45(@(t, y) S_plus(t, y, alpha), linspace(tbeg_tmp, T1, 100*T1), y0, options);
            end
            data = [data;[y(:, 1) y(:, 2)]];
            if(abs(y(end, 4)) < eps)
                % plot(y(end, 1), y(end, 2), '.r', MarkerSize=10)
            end
            if(abs(t(end) - T1) < eps)
                % plot(y(end, 1), y(end, 2), '.b', MarkerSize=10)
            end

            if sg == 1
                %plot(y(:, 1), y(:, 2), 'Color', [160,32,58]/255);
                %траектрии S+ бордовые
            elseif sg == -1
                %plot(y(:, 1), y(:, 2), 'Color', [95,223,197]/255); 
                %траектрии S- зелёные
            end
            
            if abs(t(end) - T1) < eps
                break
            end
            tbeg_tmp = t(end);
            y0 = transpose(y(end, :));
            if sg == 1
                curve2((lap - 1) * (N - 1) + i - 1, :) = y(end, 1:2);
            end
            if sg == -1
                curve1((lap - 1) * (N - 1) + i - 1, :) = y(end, 1:2);
            end
            % plot(y(end, 1), y(end, 2), '.b'); %точки переключения - синие 
            
            sg = -sg;
            lap = lap + 1;
        end
        Xfin = [Xfin; y(end, 1:2)];
    end
    
end

function [value, isterminal, direction] = events_cross_psi1(t, x)
    value = [x(4)];
    isterminal = 1;
    direction = 0;
end