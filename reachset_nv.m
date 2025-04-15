function [X, Y, X1, Y1, X2, T2] = reachset_nv(alpha, T1)

    eps = 10^-5;

    hold on
    grid on

    x0 = [0; 0; 1; 0];
    t0 = 0;
    t2 = 0;
    sg = 1;
    Xm = [];
    Ym = [];
    while 1
        if sg == 1
            options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_minus);
            [t, y] = ode45(@(t, y) S_minus(t, y, alpha), linspace(t0, T1, 500), x0, options);
        else
            options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_plus);
            [t, y] = ode45(@(t, y) S_plus(t, y, alpha), linspace(t0, T1, 500), x0, options);
        end
        Xm = [Xm; y(:, 1)];
        Ym = [Ym; y(:, 2)];
        %p = plot(y(:, 1), y(:, 2), 'r');
        %p.LineWidth = 2;
        if abs(t(end) - T1) < eps
            break;
        end
        if t2 == 0
            t2 = t(end);
        end
        %plot(y(end, 1), 0, '*')
        x0 = y(end, :);
        x0(2) = 0;
        x0(4) = 0;
        t0 = t(end);
        sg = -sg;
    end

    x0 = [0; 0; -1; 0];
    t0 = 0;
    sg = -1;
    t1 = 0;
    Xp = [];
    Yp = [];
    while 1
        if sg == 1
            options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_minus);
            [t, y] = ode45(@(t, y) S_minus(t, y, alpha), linspace(t0, T1, 500), x0, options);
        else
            options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_plus);
            [t, y] = ode45(@(t, y) S_plus(t, y, alpha), linspace(t0, T1, 500), x0, options);
        end
        if t1 == 0
            t1 = t(end);
        end
        %p = plot(y(:, 1), y(:, 2), 'b');
        %p.LineWidth = 2;
        Xp = [Xp; y(:, 1)];
        Yp = [Yp; y(:, 2)];
        %plot(0, 0, '*')
        if abs(t(end) - T1) < eps
            break;
        end
        %plot(y(end, 1), 0, '*')
        x0 = y(end, :);
        x0(2) = 0;
        x0(4) = 0;
        t0 = t(end);
        sg = -sg;
    end
    X1 = [rot90(rot90(Xm)); Xp];
    Y1 = [rot90(rot90(Ym)); Yp];
    if abs(t1) < eps
        t1 = T1;
    end
    if abs(t2) < eps
        t2 = T1;
    end


    N = 200;
    tspan1 = linspace(0, t1, N+1);
    x0 = [0; 0];
    [tbeg1, ybeg1] = ode45(@(t, y) x_plus(t, y, alpha), tspan1(1:(end-1)), x0);

    Xend1 = [];

    for i = 2:N
        x0 = [ybeg1(i, 1); ybeg1(i, 2); 1; 0];
        tbeg = tbeg1(i);
        sg = -1;
        while 1
            if sg == -1
                options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_psi);
                [t, y] = ode45(@(t, y) S_minus(t, y, alpha), [tbeg, T1], x0, options);
            else
                options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_psi);
                [t, y] = ode45(@(t, y) S_plus(t, y, alpha), [tbeg, T1], x0, options);
            end
            %plot(y(:, 1), y(:, 2));
            if abs(t(end) - T1) < eps
                break
            end
            tbeg = t(end);
            x0 = transpose(y(end, :));
            sg = -sg;
        end
        Xend1 = [Xend1; y(end, 1:2)];
    end
    %p = plot(Xend1(:, 1), Xend1(:, 2), 'g')
    %p.LineWidth = 2;
    Xend2 = [];
    tspan2 = linspace(0, t2, N+1);
    x0 = [0; 0];
    [tbeg2, ybeg2] = ode45(@(t, y) x_minus(t, y, alpha), tspan2(1:(end-1)), x0);

    for i = 1:N
        x0 = [ybeg2(i, 1); ybeg2(i, 2); 1; 0];
        tbeg = tbeg2(i);
        sg = 1;
        while 1
            if sg == -1
                options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_psi);
                [t, y] = ode45(@(t, y) S_minus(t, y, alpha), [tbeg, T1], x0, options);
            else
                options = odeset('RelTol' , 1e-2, 'AbsTol' , 1e-5, 'Events', @events_psi);
                [t, y] = ode45(@(t, y) S_plus(t, y, alpha), [tbeg, T1], x0, options);
            end
            %plot(y(:, 1), y(:, 2));
            if abs(t(end) - T1) < eps
                break
            end
            tbeg = t(end);
            x0 = transpose(y(end, :));
            sg = -sg;
        end
        Xend2 = [Xend2; y(end, 1:2)];
    end
    %disp(Xend2)

   %p = plot(Xend2(:, 1), Xend2(:, 2), 'y');
   %p.LineWidth = 2;
    
    Xh = [Xend1; Xend2];
    
    X = Xh(:, 1);
    Y = Xh(:, 2);
    D = [X Y];
    
    [X, Y] = selfIntersections(X, Y);

    X = [X1(1); X; X1(1)];
    Y = [Y1(1); Y; Y1(1)];
    
    %{
    u = alpha;
    f = @(x) -beta*x.*sin(x.^3) + u;
    
    x1 = min(X);
    x2 = 0;
    N = 1000;
    x0 = linspace(x1,x2,N);
    z1 = fzero(f, x0(1));
    j = 1;
    for i = 2:N
        z0 = fzero(f, x0(i));
        if abs(z0 - z1) > eps & z0 < 0 & z0 > x1
            j = j + 1;
            z1(j) = z0;
        end
    end
    z1 = fliplr(sort(z1));
    
    u = -alpha;
    f = @(x) - beta*x.*sin(x.^3) + u;
    
    x1 = 0;
    x2 = max(X);
    N = 500;
    x0 = fliplr(linspace(x1,x2,N));
    z2 = fzero(f, x0(1));
    j = 1;
    for i = 2:N
        z0 = fzero(f, x0(i));
        if abs(z0 - z2) > eps & z0 > 0 & z0 < x2
            j = j + 1;
            z2(j) = z0;
        end
    end
    z2 = fliplr(sort(z2));
    
    
    X2 = [z1 z2];
    T2 = [length(z1) length(z2)];
    %}
    X2 = [];
    T2 = [];
    
end

function [value, isterminal, direction] = events_psi(t, y)
    value = [y(4)];
    isterminal = 1;
    direction = 0;
end

function [value, isterminal, direction] = events_minus(t, y)
    value = [y(2)];
    isterminal = 1;
    direction = 0;
end

function [value, isterminal, direction] = events_plus(t, y)
    value = [y(2)];
    isterminal = 1;
    direction = 0;
end