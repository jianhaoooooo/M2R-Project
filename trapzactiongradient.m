function [A G] = trapzactiongrade(x, A_sym ,G_sym, N, ck)
    disp('Optimizing...')
    %syms(sym('u', [1 N]));
    %u = (sym('u',[1 N]));
    %syms(sym('v', [1 N]));
    %v = (sym('v',[1 N]));

    %ck = [u v];  
    
    A = double(subs(A_sym, ck, x));
    G = double(subs(G_sym, ck, x));
end