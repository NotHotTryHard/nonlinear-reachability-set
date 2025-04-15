%%
% main
clear, clc;

t_0 = 0;
T = 2;
alpha = 0.7;

hold on;
[X, Y, X1, Y1, X2, T2] = reachset(alpha, T);

plot(X, Y, 'Color', 'b', 'LineWidth', 2);
plot(X1, Y1, '.', 'Color', 'r', 'LineWidth', 2);

%границы для красивого вывода
left_x = min(X) - 1;
right_x = max(X) + 1;
left_y = min(Y) - 1;
right_y = max(Y) + 1;

%отображение всех точек
plot(X2(1, T2 == 1), X2(2, T2 == 1), '.', 'Color', 'g', 'MarkerSize', 10);
plot(X2(1, T2 == 2), X2(2, T2 == 2), '.', 'Color', 'black', 'MarkerSize', 10);

%отображение соответствующих точек
%First = X2(:, T2 == 1 & X2(1, :) < 0);
%Second = X2(:, T2 == 2 & X2(1, :) > 0);
%plot(First(1, :), First(2, :), '.', 'Color', 'g', 'MarkerSize', 8);
%plot(Second(1, :), Second(2, :), '.', 'Color', 'black', 'MarkerSize', 8);

% настройка выходного окна
xlim([left_x, right_x]);
ylim([left_y, right_y]);
xlabel('x_1');
ylabel('x_2');

hold off; 
%%
clear; clc

%Вывод анимации в файл и на экран

reachsetdyn(0.2, 1, 5, 5, [])
%reachsetdyn(1, 0.4, 2, 5, [])
%reachsetdyn(1, 0.4, 2, 5, 'animation.avi')
%reachsetdyn(0.5, 1, 1, 1, [])

%%

%функция события, когда x_2 обнулилась
function [value, isterminal, direction] = VEvent(t, x)
    value = x(2);  
    isterminal = 1;   
    direction = 0;   
end

%построение границы множества достижимости
function [X, Y, X1, Y1, X2, T2] = reachset(alpha, T)

% инициализация
    X = [];
    Y = [];
    X1 = [];
    Y1 = [];
    X2 = [];
    T2 = [];

% построение куска решения S_{+} до первого обнуления x_2    
    dxdt_1 = @(t, x) [x(2);
                    -3*x(1)^3 - 2*sin(x(1)) - x(2)^2 + 0.7];
    x0 = [0, 0];
    %tspan = [0, T];
    tspan = [0 : 0.025 : T];    
    [t, x] = ode45(dxdt_1, tspan, x0, odeset('Events', @VEvent));   

% итерационный процесс построения граничных точек
    t_switch = t(2 : end);
    q = 1;
    x_0 = x(2 : end, 1);
    y_0 = x(2 : end, 2);
    X1 = [0; x_0];
    Y1 = [0; y_0];
    
    while any(t_switch - T < 0)
            
        for i = 1 : size(t_switch, 1)
            if (mod(q, 2) == 1) && (t_switch(i) < T)
                [t_switch(i), x_0(i), y_0(i)] = negativeSolve(t_switch(i), x_0(i), y_0(i), T, alpha);
            else if (t_switch(i) < T)
                    [t_switch(i), x_0(i), y_0(i)] = positiveSolve(t_switch(i), x_0(i), y_0(i), T, alpha);
                end
            end
        end

        X1 = [X1; x_0(t_switch - T < 0)];
        Y1 = [Y1; y_0(t_switch - T < 0)];
        q = q + 1;
            
    end
    X = [X; x_0];
    Y = [Y; y_0];        

% аналогично для S_{-}
    dxdt_2 = @(t, x) [x(2);
            -x(1)^3 - 2*sin(x(1)) - x(2)^2 - 0.7];
    x0 = [0, 0];
    %tspan = [0, T];
    tspan = [0 : 0.025 : T];
    [t, x] = ode45(dxdt_2, tspan, x0, odeset('Events', @VEvent));   
    
    t_switch = t(2 : end);
    q = 0;
    x_0 = x(2 : end, 1);
    y_0 = x(2 : end, 2);
    X11 = [0; x_0];
    Y11 = [0; y_0];

    while any(t_switch - T < 0)

        for i = 1 : size(t_switch, 1)
            if (mod(q, 2) == 1) && (t_switch(i) < T)
                [t_switch(i), x_0(i), y_0(i)] = negativeSolve(t_switch(i), x_0(i), y_0(i), T, alpha);
            else if (t_switch(i) < T)
                    [t_switch(i), x_0(i), y_0(i)] = positiveSolve(t_switch(i), x_0(i), y_0(i), T, alpha);
                end
            end
        end

        X11 = [X11; x_0(t_switch - T < 0)];
        Y11 = [Y11; y_0(t_switch - T < 0)];
        q = q + 1;

    end
    X = [X; x_0];
    Y = [Y; y_0];
    X1 = [rot90(rot90(X11)); X1];
    Y1 = [rot90(rot90(Y11)); Y1];
    X = [X; X(1)];
    Y = [Y; Y(1)];
    
    %X11 = X1;
    %Y11 = Y1;

% удаляем самопересечения
    [X, Y] = beauty(X, Y);

    
    %[t, x] = ode45(dxdt_1, tspan, x0);
    %X1 = [0; x(:, 1)];
    %Y1 = [0; x(:, 2)];

    %[t, x] = ode45(dxdt_2, tspan, x0);
    %X1 = [rot90(rot90(x(:, 1))); X1];
    %Y1 = [rot90(rot90(x(:, 2))); Y1];  


% ищем особые точки, те, которые потребуется вывести
    left = min(X) - alpha - 1;
    right = max(X) + alpha + 1;
    left = min(left, -right);
    right = max(right, -left);

    s0 = [left : 0.1 : right];
    eps = 0.001;
    f = @(x) x .* cos(x.^2) - alpha;
    for i = 1 : size(s0, 2)
        options = optimoptions('fsolve','Display','off');
        a = fsolve(f, s0(i), options);
        if abs(a * cos(a^2) - alpha) < eps
            X2(end + 1) = a;
        end
    end
    X2 = round(X2, 4);
    X2 = unique(X2);
    T2 = ones(1, size(X2, 2));



    tmp = -X2;
    X2 = [X2, tmp];
    T2 = [T2, ones(1, size(tmp, 2)) + 1];

    X2 = [X2; zeros(1, size(X2, 2))];
    
    %X1 = X11;
    %Y1 = Y11;
end

%удаление самопересечений
function [X, Y] = beauty(X, Y)

    i = 1;
    
    while i <= size(X, 1) - 1
        
        A = [X(i), 1; X(i + 1), 1];
        B = [Y(i); Y(i+1)];
        C = linsolve(A, B);    
        j = 1;
        while j <= size(X, 1) - 1        
            if (i ~= j) && (i ~= j - 1) && (i ~= j + 1)
                
                res_1 = C(1)*X(j) + C(2) - Y(j);
                res_2 = C(1)*X(j+1) + C(2) - Y(j+1);
                
                if (res_1 * res_2 < 0)
                    
                    A = [X(j), 1; X(j + 1), 1];
                    B = [Y(j); Y(j+1)];
                    D = linsolve(A, B);                    
                    res_3 = D(1)*X(i) + D(2) - Y(i);
                    res_4 = D(1)*X(i+1) + D(2) - Y(i+1);
                    
                    if (res_3 * res_4 < 0)
                        A = [C(1), -1; D(1), -1];
                        B = [-C(2); -D(2)];
                        E = linsolve(A, B);
                        
                        X(i + 1 : j) = [];
                        X = [X(1 : i); E(1); X(i+1 : end)];
                        Y(i + 1 : j) = [];
                        Y = [Y(1 : i); E(2); Y(i+1 : end)];
                    end
                end
            end
            j = j + 1;
        end
        i = i + 1;
    end
end

%функция события, когда пси_2 обнуляется
function [value, isterminal, direction] = PsiEvent(t, x)
    value = x(4);  
    isterminal = 1;   
    direction = 0;   
end

%решение системы S_{+}
function [t_switch, x_1, y_1] = positiveSolve(t_switch, x_0, y_0, T, alpha)
    dxdt = @(t, x) [x(2);
            -3*x(1)^3 - 2*sin(x(1)) - x(2)^2 + 0.7;
            x(4)*(-9*x(1)^2 - 2*cos(x(1)));
            x(3) -2*x(2)*x(4)];     
           
    x0 = [x_0, y_0, -1, 0];
    tspan = linspace(t_switch, T, T*100);
    
    [t, x] = ode45(dxdt, tspan, x0, odeset('Events', @PsiEvent));
    x_1 = x(end, 1);
    y_1 = x(end, 2);
    t_switch = t(end);
    plot(x(:, 1), x(:, 2), 'Color', 'green');
end

%решение системы S_{-}
function [t_switch, x_1, y_1] = negativeSolve(t_switch, x_0, y_0, T, alpha)
    dxdt = @(t, x) [x(2);
            -3*x(1)^3 - 2*sin(x(1)) - x(2)^2 - 0.7;
            x(4)*(-9*x(1)^2 - 2*cos(x(1)));
            x(3) -2*x(2)*x(4)]; 
           
    x0 = [x_0, y_0, 1, 0];
    tspan = linspace(t_switch, T, T*100);
    [t, x] = ode45(dxdt, tspan, x0, odeset('Events', @PsiEvent));
    x_1 = x(end, 1);
    y_1 = x(end, 2);
    t_switch = t(end);
    plot(x(:, 1), x(:, 2), 'Color', 'black');
end

%анимация в файл и вывод на экран
function reachsetdyn(alpha, t1, t2, N, filename)

    t = linspace(t1, t2, N);
    
    if filename

        hold on;
        xlabel('x_1');
        ylabel('x_2');
        legend('Location', 'southeast');
        str = strcat('\alpha = ', ' ', num2str(alpha));        
        title(str);
        for i = 1 : N 
            [X, Y] = reachset(alpha, t(i));
            str = strcat('t = ', num2str(t(i)));
            plot(X, Y, 'LineWidth', 2, 'DisplayName', str);
            frame(i) = getframe(gcf)           
        end
        hold off;

        v = VideoWriter(filename);
        v.FrameRate = 1/3;
        open(v)
        writeVideo(v,frame)
        close(v)
        
    else
        
        hold on;
        xlabel('x_1');
        ylabel('x_2');
        legend('Location', 'southeast');
        str = strcat('\alpha = ', ' ', num2str(alpha));        
        title(str);     
        for i = 1 : N 
            [X, Y] = reachset(alpha, t(i));
            str = strcat('t = ', num2str(t(i)));
            plot(X, Y, 'LineWidth', 2, 'DisplayName', str);
            pause(3);
        end
        hold off;
        
    end
    
end
