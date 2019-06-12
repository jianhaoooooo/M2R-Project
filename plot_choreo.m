function graph = plot_choreo(n, N, c)
    syms t
    k = -(N-1)/2:(N-1)/2;
    tp = linspace(0, 2*pi);
    complex_c = c(1:N) + 1i*c(N+1:2*N); 
    complex_c = reshape(complex_c,[1, N]);
    choreo = sum(complex_c.*exp(1i.*k.*t));
    choreo_graph = subs(choreo, t, tp);
    
    dots = linspace(0, 2*pi*(n-1)/n, n);
    eq_points = double(subs(choreo, t, dots(1:3)));
    
    plot(choreo_graph, '-', 'LineWidth', 1.5)
    hold on
    plot(eq_points, '.', 'MarkerSize', 30)
    hold off
end