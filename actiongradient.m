function [A, G] = actiongradient(x, A_sym ,G_sym, ck)
    disp('Optimizing...')
    A = double(subs(A_sym, ck, x));
    
    if nargout >1 
        G = double(subs(G_sym, ck, x));
end