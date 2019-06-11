function graph = plot_choreo(N, c)
    syms t
    tp = linspace(0, 2*pi);
    k = -(N-1)/2:(N-1)/2;
    complex_c = c(1:N) + 1i*c(N+1:2*N); 
    complex_c = reshape(complex_c,[1, N]);
    choreo = sum(complex_c.*exp(1i.*k.*t));
    choreo = subs(choreo, t, tp);
    plot(choreo)
end