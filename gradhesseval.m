% c needs to be a 2M vector [u1 u2 .. uM v1 v2 ... vM]

function [G, H] = gradhesseval(c, n, M, w) 

syms(sym('u', [1 M]))
u = (sym('u',[1 M]));
syms(sym('v', [1 M]))
v = (sym('v',[1 M]));

% grad 
k = -((M-1)/2) : ((M-1)/2);
% grad kinetic part
Ak_uk = 2*sym('pi')*n.*((k+w).^2).*u; 
Ak_vk = 2*sym('pi')*n.*((k+w).^2).*v;
grad_Ak = [Ak_uk Ak_vk];
grad_Ak = subs(grad_Ak, [u v], c); % need to make sure what c actually is, is it [u1 ... uM v1... vM]?

%grad potential part


grad_Ap = sym(zeros(2*M, 1));
for j=1:n-1 % summation over j
    k = -((M-1)/2) : ((M-1)/2);
    syms t real
    a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
    b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);
    f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
    power_f_j = f_j^(1.5);
    
    %grad_f_j
    l = -((M-1)/2) : ((M-1)/2);
    grad_f_j = sym(zeros(2*M, 1));

    for k=-((M-1)/2) : ((M-1)/2)
        syms t real
        a_lj = (1-cos(2*sym('pi')*j.*l./n)).*cos(l.*t) + sin(2*sym('pi').*j.*l./n).*sin(l.*t);
        b_lj = (-1+cos(2*sym('pi')*j.*l./n)).*sin(l.*t) + sin(2*sym('pi').*j.*l./n).*cos(l.*t);
        a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
        b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);

        tmp = 2 * sum((a_lj.*a_kj + b_lj.*b_kj).*u + (b_lj.*a_kj - a_lj.*b_kj).*v);
        grad_f_j(k+(M+1)/2) = tmp;
    end

    for k=-((M-1)/2) : ((M-1)/2)
        syms t real
        a_lj = (1-cos(2*sym('pi')*j.*l./n)).*cos(l.*t) + sin(2*sym('pi').*j.*l./n).*sin(l.*t);
        b_lj = (-1+cos(2*sym('pi')*j.*l./n)).*sin(l.*t) + sin(2*sym('pi').*j.*l./n).*cos(l.*t);
        a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
        b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);

        tmp = 2 * sum((a_lj.*b_kj - b_lj.*a_kj).*u + (b_lj.*b_kj - a_lj.*a_kj).*v);
        grad_f_j(k+(3*M+1)/2) = tmp;
    end
    
    expression_inside = grad_f_j./power_f_j; % transposed grad fj
    expression_inside = subs(expression_inside, [u v], c); % need to make sure what c actually is, is it [u1 ... uM v1... vM]?
    integral_part = vpaintegral(expression_inside, t, 0, 2*pi); % int takes extremely long, vpaintegral works quickly but it's numerical integration
    final_expr_j = -(n/4).*integral_part;
    grad_Ap = grad_Ap + final_expr_j; %summation over j via loop
    
end

G = grad_Ak' + grad_Ap;

% hessaval
H = sym(zeros(2*M, 2*M));

% kinetic part
for k=-((M-1)/2) : ((M-1)/2)
    H(k+(M+1)/2, k+(M+1)/2) = 2*pi*n*(k^2);
    H(k+(3*M+1)/2, k+(3*M+1)/2) = 2*pi*n*(k^2);
end

% potential part
Au_ul_uk = sym(zeros(M, M));
Au_vl_vk = sym(zeros(M, M));
Au_ul_vk = sym(zeros(M, M));

for j=1: n-1
    k = -((M-1)/2) : ((M-1)/2);
    l = -((M-1)/2) : ((M-1)/2);
    syms t real
    a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
    b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);
    f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
    
    % for Au_ul_uk
    fj_ul = G(1:M);
    product_fj_ul = fj_ul*fj_ul';
    
    fj_ul_uk = 2.*(a_kj' * a_kj + b_kj' * b_kj);
    
    internal_int1 = (fj_ul_uk .* f_j - (3/2) .* product_fj_ul)./(f_j)^(2.5);
    internal_int1 = subs(internal_int1, [u v], c);
    
    tmp = vpaintegral(internal_int1, t, 0, 2*pi);
    Au_ul_uk = Au_ul_uk + -n/4 .* tmp;
    
    % similarly for Au_vl_vk
    fj_vl = G(M+1:end);
    product_fj_vl = fj_vl*fj_vl';
    fj_vl_vk = fj_ul_uk;
    
    internal_int2 = (fj_vl_vk .* f_j - (3/2) .* product_fj_vl)./(f_j)^(2.5);
    internal_int2 = subs(internal_int2, [u v], c);
    
    tmp2 = vpaintegral(internal_int2, t, 0, 2*pi);
    Au_vl_vk = Au_vl_vk + -n/4 .* tmp2;
    
    % for Au_ul_vk
    fj_ul = G(1:M);
    fj_vl = G(M+1:end);
    product_fj_ul_vk = fj_ul*fj_vl';
    fj_ul_vk = 2.*(a_kj' * b_kj - b_kj' * a_kj);
    
    internal_int3 = (fj_ul_vk .* f_j - (3/2) .* product_fj_ul_vk)./(f_j)^(2.5);
    internal_int3 = subs(internal_int3, [u v], c);
    
    tmp3 = vpaintegral(internal_int3, t, 0, 2*pi);
    Au_ul_vk = Au_ul_vk + -n/4 .* tmp3;
end

H(1:M, 1:M) = H(1:M, 1:M) + Au_ul_uk;
H(M+1:end, M+1:end) = H(M+1:end, M+1:end) + Au_vl_vk;
H(M+1:end, 1:M) = H(M+1:end, 1:M) + Au_ul_vk;
H(1:M, M+1:end) = H(M+1:end, 1:M)';

% setting all the correct zeros (to make sure)
H((M+1)/2 , :) = 0;
H(: , (M+1)/2) = 0;
H(: , ((3*M-1)/2)) = 0;
H(((3*M-1)/2), :) = 0;

end