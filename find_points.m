range = 0.2;
mesh = linspace(-range, range, 100);
alpha = 0.1;
f1 = @(x) 3*x.*x.*x +2 * sin(x) - alpha;
f2 = @(x) 3*x.*x.*x +2 * sin(x) + alpha;

pt1 = fzero(f1, 0)
pt2 = fzero(f2, 0)

p1 = plot(mesh, f1(mesh), 'color', [160,32,58]/320,LineWidth=1.2);
hold on;
grid on;
xlabel('x')
ylabel('f(x)')

ax = gca;
p2 = plot(mesh, f2(mesh), 'color', [95,223,197]/320, LineWidth=1.2);
plot(linspace(-range, range, 20), repmat(0, 1, 20),'k', LineWidth=1.2);
plot(pt1, 0,'pentagram', 'MarkerSize',10, 'MarkerEdgeColor',[160,32,58]/255,'MarkerFaceColor', [160,32,58]/255) % особая точка S-

plot(pt2, 0,'pentagram','MarkerSize',10,'MarkerEdgeColor',[95,223,197]/255,'MarkerFaceColor', [95,223,197]/255) % особая точка S+



legend([p1 p2], {'(S_+): y = 3x^3+2sin(x)-'+string(alpha), '(S_-): y = 3x^3+2sin(x)+'+string(alpha)})
exportgraphics(ax,'plot.eps')