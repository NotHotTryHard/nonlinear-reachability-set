function [XX, YY] = selfIntersections(X, Y)
    i = 1;
    eps = 10^-10;
    while i < length(X)
        %disp('i = ')
        %disp(i)
        line1_dot1 = [X(i) Y(i)];
        line1_dot2 = [X(i+1) Y(i+1)];
        if abs(line1_dot2(1)-line1_dot1(1)) > eps
            w1 = (line1_dot2(2)-line1_dot1(2))/(line1_dot2(1)-line1_dot1(1));
            A1 = w1;
            B1 = -1;
            C1 = -line1_dot1(1)*w1+line1_dot1(2);
        else
            A1 = 1;
            B1 = 0;
            C1 = -line1_dot1(1);
        end
        %line1 = @(x) A1*x(1) + B1*x(2) + C1;
        %disp([A1 B1 C1])
        j = 1;
        while j < length(Y)
            %disp('j = ')
            %disp(j)
            if abs(i-j) <= 1
                j = j + 1;
                continue
            end
            line2_dot1 = [X(j) Y(j)];
            line2_dot2 = [X(j+1) Y(j+1)];
            %disp("ajpfe")
            %disp(A1*line2_dot1(1) + B1*line2_dot1(2) + C1)
            %disp(A1*line2_dot2(1) + B1*line2_dot2(2) + C1)
            if (A1*line2_dot1(1) + B1*line2_dot1(2) + C1)*(A1*line2_dot2(1) + B1*line2_dot2(2) + C1) > 0
                j = j + 1;
                continue;
            end
            %disp([A1 B1 C1])
%             disp(A1*line2_dot1(1) + B1*line2_dot1(2) + C1)
%             disp(A1*line2_dot2(1) + B1*line2_dot2(2) + C1)
%             disp([i j])
            %disp(line1_dot1)
            %disp(line1_dot2)
            %disp(line2_dot1)
            %disp(line2_dot2)
            
            if abs(line2_dot2(1)-line2_dot1(1)) > eps
                w2 = (line2_dot2(2)-line2_dot1(2))/(line2_dot2(1)-line2_dot1(1));
                A2 = w2;
                B2 = -1;
                C2 = -line2_dot1(1)*w2+line2_dot1(2);
                %disp([A2 B2 C2])
                y = (A1*C2-A2*C1)/(B1*A2-A1*B2);
                x = -B1*y/A1 - C1/A1;
            else
                x = line2_dot2(1);
                y = C1 + A1*x;
            end
            if (x >= min(line1_dot1(1), line1_dot2(1))) & (x <= max(line1_dot1(1), line1_dot2(1)))
                X(i+1) = x;
                Y(i+1) = y;
                X((i+2):j) = [];
                Y((i+2):j) = [];
                j = i + 2;
            else
                j = j + 1;
            end
        end
        i = i + 1;
    end
    XX = X;
    YY = Y;
end
            
            
            