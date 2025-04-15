function [X, Y, X1, Y1, X2, T2] = reachset(alpha, T1)

    eps = 10^-5;

    hold on
    grid on

    [t1, X_u0_pos, Y_u0_pos] = find_cross_x1(alpha,[0 0 1 0]',T1, 1);
    [t2, X_u0_neg, Y_u0_neg] = find_cross_x1(alpha,[0 0 -1 0]',T1, -1);
    


    
    
    if abs(t1) < eps
        t1 = T1;
    end
    if abs(t2) < eps
        t2 = T1;
    end


    N = 400;
    curve1 = [0, 0];
    curve2 = [0, 0];
    [Xend1, data1, curve1, curve2] = find_cross_psi1(alpha, @x_plus, t1, T1, -1, N, curve1, curve2);
    [Xend2, data2, curve1, curve2] = find_cross_psi1(alpha, @x_minus, t2, T1, 1, N, curve1, curve2);
    % curve1
    % curve2 
    X1 = [rot90(rot90(X_u0_neg)); X_u0_pos];
    Y1 = [rot90(rot90(Y_u0_neg)); Y_u0_pos];
    X1 = [rot90(rot90(curve1(:,1)));rot90(rot90(X_u0_neg)); X_u0_pos; curve2(:,1)];
    Y1 = [rot90(rot90(curve1(:,2)));rot90(rot90(Y_u0_neg)); Y_u0_pos; curve2(:,2)];
    % [X1, Y1]
    if curve1 == [0 0]
        X1 = X1(2:end,:);
        Y1 = Y1(2:end,:);
    end
    if curve2 == [0 0]
        X1 = X1(1:end-1,:);
        Y1 = Y1(1:end-1,:);
    end
    % [X1, Y1]
    %{
    data = [data1; data2];
    qwe = data(:, 1);
    rty = data(:, 2);
    k = boundary(qwe,rty, 1);
    a = qwe(k);
    b = rty(k);
    k = boundary(a, b);
    a = a(k);
    b = b(k);
    %}

    Xh = [Xend1; Xend2];
    
    X = Xh(:, 1);
    Y = Xh(:, 2);
    D = [X Y];
    
    [X, Y] = selfIntersections(X, Y);

    X = [X; X(1)];
    Y = [Y; Y(1)];

    u = alpha;
    f = @(x) 3*x.^3 + 2*sin(x) - u;
    g = @(x) 3*x.^3 + 2*sin(x) + u;
    p1 = fzero(f, 0);
    p2 = fzero(g, 0);
    
    X2 = [p1 p2]; % первая - S+, вторая - S-
    T2 = [1 2]; 

    % X = a;
    % Y = b;
    % k = boundary(X, Y);
    % X = X(k);
    % Y = Y(k);
    
end
