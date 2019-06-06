% n -> number of bodies
% N -> number of BFGS Fourier coefficients
% M -> number of Newton Fourier coefficients
% if absolute choreography simply set w to 0

% choices of n, N (must be odd!), M, w 
n = 3; N = 5; M = 145; w=0;
q0 = chebfun(@(t)cos(t)+1i*sin(2*t), [0 2*pi], N,'trig'); 
c0 = trigcoeffs(q0);

%%
% kinetic part of action
syms(sym('u', [1 N]))
u = (sym('u',[1 N]));
syms(sym('v', [1 N]))
v = (sym('v',[1 N]));

k = -((N-1)/2) : ((N-1)/2);
Ak_sym = sym('pi')*n * sum(k.*(u.^2 + v.^2));

x=transpose([u,v]);

% verifying that this is indeed agrees with the explicit representation in
% the book
gradf = jacobian(Ak_sym,x).';

% use this to convert symbolic function to anonymous function Ak = matlabFunction(Ak_sym, 'vars', {[u v]}); 
options = optimoptions('fminunc');
options.GradObj = 'off';
options.Algorithm = 'quasi-newton';
options.HessUpdate = 'bfgs';
Ak = matlabFunction(Ak_sym, 'vars', {[u v]}); 
c = fminunc(Ak,[real(c0);imag(c0)]',options);

% potential part of action

Au_sym = 0;

for j= 1:n-1
    syms t real
    a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
    b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);
    f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
    tmp = (f_j)^-0.5;
    
    % this part won't evaluate properly
    integral = int(tmp, t, 0, 2*pi);
    
    % summation over j
    Au_sym = Au_sym + n/2 * integral;
end 

% the parts here won't run properly unless the integral is evaluated
% in order to have a Au_sym without the variable t


% combining kinetic and potential part together 
A_sym = Ak_sym + Au_sym;
A = matlabFunction(A_sym, 'vars', {[u v]});


options = optimoptions('fminunc');
options.GradObj = 'on';
options.Algorithm = 'quasi-newton';
options.HessUpdate = 'bfgs';
c = fminunc(@(x) actiongradeval(x, n),[real(c0);imag(c0)]',options);


%%
% Newton method:
c = [zeros((M-N)/2,1);c(1:N);zeros(M-N,1);c(N+1:end);zeros((M-N)/2,1)]; mid = 1 + floor(M/2);
[G, H] = gradhesseval(c,n); [L, D] = ldl(H);

for k = 1:2 % specify the number of iterations for Newton Method
s = L'\(D\(L\(-G)));
cnew = c + [s(1:mid-1);0;s(mid:M+mid-2);0;s(M+mid-1:end)]; c = cnew; G = gradhesseval(c,n);
end

%%
% Reconstruct solution:
c = c(1:M) + 1i*c(M+1:2*M);
q = chebfun(c,[0 2*pi],'coeffs','trig');



