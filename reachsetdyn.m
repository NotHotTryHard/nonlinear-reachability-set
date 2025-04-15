function a = reachsetdyn(alpha, t1, t2, N, filename)
    eps = 10^-6;
    if abs(t1-t2) < eps
        t = t1;
        N = 0;
    else
        t = linspace(t1, t2, N+1);
    end
    nFrames = N+1;
    mov(1:nFrames) = struct('cdata', [], 'colormap', []);
    
    [X, Y, X1, Y1] = reachset(alpha, t2);
    limx = [min(X)-0.1 max(X)+0.1];
    limy = [min(Y)-0.1 max(Y)+0.1];
    
    for i = 1:(N+1)
        [X, Y, X1, Y1] = reachset(alpha, t(i));
        hold on
        p1= plot(X, Y, 'k');
        p1.LineWidth = 1;
        p2 = plot(X1, Y1, 'Color', [255,99,71]/255);
        p2.LineWidth = 3;
        
        hold off
        xlabel('x1')
        ylabel('x2')
        xlim(limx)
        ylim(limy)
        legend([p1 p2], {'Граница области достижимости', 'Кривая переключений'})
        mov(i) = getframe();
    end
    %movie(mov,1,1)
    if filename ~= ""
        wr = VideoWriter(filename);
        wr.FrameRate = 2;
        open(wr)
        for i = 1:nFrames
            writeVideo(wr, mov(i))
        end
        close(wr);
    end
    a = 2;
end