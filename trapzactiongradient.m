function [A G] = trapzactiongradient(x, A_sym ,G_sym, N, ck)
    disp('Optimizing...')
    A = double(subs(A_sym, ck, x));
    G = double(subs(G_sym, ck, x));
end