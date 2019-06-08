% n -> number of bodies
% N -> number of BFGS Fourier coefficients
% M -> number of Newton Fourier coefficients
% if absolute choreography simply set w to 0

% choices of n, N (must be odd!), M, w 
n = 3; N = 5; M = 145; w=0;
q0 = chebfun(@(t)cos(t)+1i*sin(2*t), [0 2*pi], N,'trig'); 
c0 = trigcoeffs(q0);

%%
options = optimoptions('fminunc');
options.GradObj = 'on';
options.Algorithm = 'quasi-newton';
options.HessUpdate = 'bfgs';

objfunction = actiongradeval(n, N, w)
c = fminunc(@(x) objfunction(x),[real(c0);imag(c0)]',options);

%%
% Newton method:
c = [zeros((M-N)/2,1);c(1:N);zeros(M-N,1);c(N+1:end);zeros((M-N)/2,1)]; mid = 1 + floor(M/2);
[G, H] = gradhesseval(c,n,M,w); [L, D] = ldl(H);

for k = 1:2 % specify the number of iterations for Newton Method
s = L'\(D\(L\(-G)));
cnew = c + [s(1:mid-1);0;s(mid:M+mid-2);0;s(M+mid-1:end)]; c = cnew; G = gradhesseval(c,n);
end

%%
% Reconstruct solution:
c = c(1:M) + 1i*c(M+1:2*M);
q = chebfun(c,[0 2*pi],'coeffs','trig');



