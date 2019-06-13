% n -> number of bodies
% N -> number of BFGS Fourier coefficients
% M -> number of Newton Fourier coefficients
% if absolute choreography simply set w to 0

% choices of n, N (must be odd!), M, w 
n = 3; N = 29; M = 39; w=0;
q0 = chebfun(@(t)cos(t)+1i*sin(2*t), [0 2*pi], N,'trig'); 
c0 = trigcoeffs(q0);

%%
options = optimoptions('fminunc');
options.GradObj = 'on';
options.Algorithm = 'quasi-newton';
options.HessUpdate = 'bfgs';

[A_sym ck G_sym] = symactiongradient(n, N, w);
disp('Obtained symbolic expressions')
objfunction = @(x) actiongradient(x, A_sym, G_sym, ck);
[c fval] = fminunc(@(x) objfunction(x),transpose([real(c0);imag(c0)]),options);
c_bfgs = c; A_bfgs = fval;

%%
% Newton method:
c = reshape(c, [2*N, 1]);
c = [zeros((M-N)/2,1);c(1:N);zeros(M-N,1);c(N+1:end);zeros((M-N)/2,1)]; 
mid = 1 + floor(M/2);
[G, H] = gradienthesseval(c,n,M,w); 
[L, D] = ldl(H);

for k = 1:2 % specify the number of iterations for Newton Method
    s = L'\(D\(L\(-G)));
    new_c = c + [s(1:mid-1);0;s(mid:M+mid-2);0;s(M+mid-1:end)]; 
    c = new_c;
    G = gradienthesseval(c,n,M,w);
end

c_newton = transpose(c);
[A_sym_newton ck] = symactiongradient(n, M, w);
A_newton = actiongradient(c_newton, A_sym_newton, 0, ck) 
%%
% Constructing choreography
% input c will be a vector of u1...uk v1... vk
plot_choreo(n, M, c_newton)



