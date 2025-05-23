clc, clear
% alpha = 1.753828;

%    t  |  25  |  25  |  25  |  25  |  23  |  16  |  12  |  10  |   8  |  5.4 |  4.5 |   4  |  3.8 |  2.5 |  2.2 |  2.5 |                             2 
% alpha | 0.01 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 | 0.95 | 1.00 |   

% alpha = 1.00;
% T1 = 2.5;

alpha = 0.54;
T1 = 4.3;
eps = 10^-4;

hold on
grid on
ax = gca;
% [X, Y, X1, Y1, X2, T2, a, b, qwe, rty] = reachset(alpha, T1);
[X, Y, X1, Y1, X2, T2] = reachset(alpha, T1);

xlim([min(X)-0.1 max(X)+0.1])
ylim([min(Y)-0.1 max(Y)+0.1])
% k = boundary(a, b);
p1 = plot(X, Y, 'k');
p1.LineWidth = 1;

p2 = plot(X1, Y1,'Color', [255,99,71]/255); % кривая переключения -
% оранжевая
p2.LineWidth = 2;
% plot(qwe, rty, '.','Color', [250,250,250]/255);
%plot(a, b)


plot(0,0, 'k.', MarkerSize=8)
 
plot(X2(1), 0,'pentagram', 'MarkerSize',8, 'MarkerEdgeColor',[160,32,58]/350,'MarkerFaceColor', [160,32,58]/255) % особая точка S+
plot(X2(2), 0,'pentagram','MarkerSize',8,'MarkerEdgeColor',[95,223,197]/350,'MarkerFaceColor', [95,223,197]/255) % особая точка S-

 
xlabel('x_1')
ylabel('x_2')
exportgraphics(ax,'plot3_5.eps')

%%
N = 12;
a = reachsetdyn(alpha, 0.1, T1, N, "video3");

%%
