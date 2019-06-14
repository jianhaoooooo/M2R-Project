function graph = plot_choreo(n, N, c)
    syms t
    k = -(N-1)/2:(N-1)/2;
    tp = linspace(0, 2*pi);
    complex_c = c(1:N) + 1i*c(N+1:2*N); 
    complex_c = reshape(complex_c,[1, N]);
    choreo = sum(complex_c.*exp(1i.*k.*t));
    choreo_graph = subs(choreo, t, tp);
    
    dots = linspace(0, 2*pi*(n-1)/n, n);
    eq_points = double(subs(choreo, t, dots(1:n)));
    
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 9])
    plot(choreo_graph, '-', 'LineWidth', 9)
    axis([-1.5 1.5 -0.75 0.75])
    set(gca,'FontSize',60)
    hold on
    plot(eq_points, '.', 'MarkerSize', 130)
    hold off
end