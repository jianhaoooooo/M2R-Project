% n -> number of bodies
% N -> number of BFGS Fourier coefficients
% M -> number of Newton Fourier coefficients
% if absolute choreography simply set w to 0

% choices of n, N (must be odd!), M, w 
n = 3; N = 15; M = 145; w=0;
q0 = chebfun(@(t)cos(t)+1i*sin(2*t), [0 2*pi], N,'trig'); 
c0 = trigcoeffs(q0);

%%
options = optimoptions('fminunc');
options.GradObj = 'on';
options.Algorithm = 'quasi-newton';
options.HessUpdate = 'bfgs';

[A_sym G_sym, c_k] = actiongradient(n, N, w);
disp('Obtained symbolic expressions')
objfunction = @(x) trapzactiongradient(x, A_sym, G_sym, N, c_k);
[c fval] = fminunc(@(x) objfunction(x),transpose([real(c0);imag(c0)]),options);
c_bfgs = c; A_bfgs = fval

%%
% Newton method:
c = transpose(c);
newc = [zeros((M-N)/2,1);c(1:N);zeros(M-N,1);c(N+1:end);zeros((M-N)/2,1)]; 
mid = 1 + floor(M/2);
newc = transpose(newc);
[G, H] = gradienthesseval(newc,n,M,w); 
[L, D] = ldl(H);

for k = 1:2 % specify the number of iterations for Newton Method
s = L'\(D\(L\(-G)));
%tmp1 = L\(-G)
%tmp2 = D\tmp1
%s = L'\tmp2
cnew = newc + [s(1:mid-1);0;s(mid:M+mid-2);0;s(M+mid-1:end)]; 
c = cnew; 
G = gradhesseval(c,n,M,w);
end

%%
% Constructing choreography
%c = c(1:M) + 1i*c(M+1:2*M);
%q = chebfun(c,[0 2*pi],'coeffs','trig');

N=55
c = c_bfgs(1:N) + 1i*c_bfgs(N+1:2*N);
t = linspace(0, 2*pi);
fnew = 0
for k=1:N
    tmp = c(k)*exp(1i*(k-((N+1)/2))*t);
    fnew = fnew + tmp
end
plot(fnew)


