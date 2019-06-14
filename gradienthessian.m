% c needs to be a 2M vector [u1 u2 .. uM v1 v2 ... vM]

function [G, H] = gradienthessian(c, n, M, w)
c = reshape(c, [1, 2*M]);
u = c(1:M);
v = c(M+1:end);

% grad 
k = -((M-1)/2) : ((M-1)/2);
% grad kinetic part
Ak_uk = 2*pi*n.*((k+w).^2).*u; 
Ak_vk = 2*pi*n.*((k+w).^2).*v;
grad_Ak = [Ak_uk Ak_vk];

%grad potential part
grad_Ap = zeros(2*M, 1);
for j=1:n-1 % summation over j
    k = -((M-1)/2) : ((M-1)/2);
    syms t real
    a_kj = (1-cos(2*pi*j.*k./n)).*cos(k.*t) + sin(2*pi.*j.*k./n).*sin(k.*t);
    b_kj = (-1+cos(2*pi*j.*k./n)).*sin(k.*t) + sin(2*pi.*j.*k./n).*cos(k.*t);
    f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
    power_f_j = f_j^(1.5);
    
    %grad_f_j
    l = -((M-1)/2) : ((M-1)/2);
    grad_f_j = sym(zeros(2*M, 1));
    a_lj = (1-cos(2*pi*j.*l./n)).*cos(l.*t) + sin(2*pi.*j.*l./n).*sin(l.*t);
    b_lj = (-1+cos(2*pi*j.*l./n)).*sin(l.*t) + sin(2*pi.*j.*l./n).*cos(l.*t);

    for k=-((M-1)/2) : ((M-1)/2)
        syms t real
        a_kj = (1-cos(2*pi*j.*k./n)).*cos(k.*t) + sin(2*pi.*j.*k./n).*sin(k.*t);
        b_kj = (-1+cos(2*pi*j.*k./n)).*sin(k.*t) + sin(2*pi.*j.*k./n).*cos(k.*t);

        tmp = 2 * sum((a_lj.*a_kj + b_lj.*b_kj).*u + (b_lj.*a_kj - a_lj.*b_kj).*v);
        grad_f_j(k+(M+1)/2) = tmp;
        
        tmp2 = 2 * sum((a_lj.*b_kj - b_lj.*a_kj).*u + (b_lj.*b_kj + a_lj.*a_kj).*v);
        grad_f_j(k+(3*M+1)/2) = tmp2;
    end

    expression_inside = grad_f_j./power_f_j;
    disp('Starting integration in gradient parts')
    integral_part = zeros(2*M, 1);
    for i=1:(2*M)
        integral_part(i) = vpaintegral(expression_inside(i), t, 0, 2*pi);
        status = [num2str(i*100/(2*M)), '% completed'];
        disp(status)
    end
    
    final_expr_j = -(n/4).*integral_part;
    grad_Ap = grad_Ap + final_expr_j; %summation over j via loop
    
end

G = transpose(grad_Ak) + grad_Ap;
G = double(G);
% removing rows and columns corresponding to u0 v0
G((M+1)/2) = []; G(((3*M+1)/2)-1) = [];
disp('Completed gradient computation')

if nargout >1 
    disp('Starting hessaval computation')
    % hessaval
    Hk = sym(zeros(2*M, 2*M));

    % kinetic part
    for k=-((M-1)/2) : ((M-1)/2)
        Hk(k+(M+1)/2, k+(M+1)/2) = 2*pi*n*(k^2);
        Hk(k+(3*M+1)/2, k+(3*M+1)/2) = 2*pi*n*(k^2);
    end

    % potential part
    Au_ul_uk = zeros(M, M);
    Au_vl_vk = zeros(M, M);
    Au_ul_vk = zeros(M, M);

    for j=1: n-1
        k = -((M-1)/2) : ((M-1)/2);
        l = -((M-1)/2) : ((M-1)/2);
        
        syms t real
        a_kj = (1-cos(2*pi*j.*k./n)).*cos(k.*t) + sin(2*pi.*j.*k./n).*sin(k.*t);
        b_kj = (-1+cos(2*pi*j.*k./n)).*sin(k.*t) + sin(2*pi.*j.*k./n).*cos(k.*t);
        f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
        
        %grad_f_j
        grad_f_j = sym(zeros(2*M, 1));
        a_lj = (1-cos(2*pi*j.*l./n)).*cos(l.*t) + sin(2*pi.*j.*l./n).*sin(l.*t);
        b_lj = (-1+cos(2*pi*j.*l./n)).*sin(l.*t) + sin(2*pi.*j.*l./n).*cos(l.*t);

        for K=-((M-1)/2) : ((M-1)/2)
            syms t real
            a_Kj = (1-cos(2*pi*j.*K./n)).*cos(K.*t) + sin(2*pi.*j.*K./n).*sin(K.*t);
            b_Kj = (-1+cos(2*pi*j.*K./n)).*sin(K.*t) + sin(2*pi.*j.*K./n).*cos(K.*t);

            tmp = 2 * sum((a_lj.*a_Kj + b_lj.*b_Kj).*u + (b_lj.*a_Kj - a_lj.*b_Kj).*v);
            grad_f_j(K+(M+1)/2) = tmp;

            tmp2 = 2 * sum((a_lj.*b_Kj - b_lj.*a_Kj).*u + (b_lj.*b_Kj + a_lj.*a_Kj).*v);
            grad_f_j(K+(3*M+1)/2) = tmp2;
        end

        % for Au_ul_uk
        fj_ul = grad_f_j(1:M);
        product_fj_ul = fj_ul*transpose(fj_ul);

        fj_ul_uk = 2.*(transpose(a_kj) * a_kj + transpose(b_kj) * b_kj);

        internal_int1 = (fj_ul_uk .* f_j - (3/2) .* product_fj_ul)./((f_j)^(2.5));
        tmp = zeros(M, M);
        
        disp('Starting integration for equation (7.12) for ul uk')
        for i=1:M
            for j=1:M
            tmp(i,j) = vpaintegral(internal_int1(i,j), t, 0, 2*pi);
            end
            status = [num2str(i*100/(M)), '% completed'];
            disp(status)
        end
        
        Au_ul_uk = Au_ul_uk - n/4 .* tmp; 
        
        % for Au_vl_vk
        fj_vl = grad_f_j(M+1:end);
        product_fj_vl = fj_vl*transpose(fj_vl);
        fj_vl_vk = fj_ul_uk;
        

        internal_int2 = (fj_vl_vk .* f_j - (3/2) .* product_fj_vl)./((f_j)^(2.5));
        tmp2 = zeros(M, M);
        
        disp('Starting integration for equation (7.12) for vl vk')
        for i=1:M
            for j=1:M
            tmp2(i,j) = vpaintegral(internal_int2(i,j), t, 0, 2*pi);
            end
            status = [num2str(i*100/(M)), '% completed'];
            disp(status)
        end
        
        Au_vl_vk = Au_vl_vk - n/4 .* tmp2; 
        
        % for Au_ul_vk
        fj_ul = grad_f_j(1:M);
        fj_vl = grad_f_j(M+1:end);
        product_fj_ul_vk = fj_ul*transpose(fj_vl);
        fj_ul_vk = 2.*(transpose(a_kj) * b_kj - transpose(b_kj) * a_kj);

        internal_int3 = (fj_ul_vk .* f_j - (3/2) .* product_fj_ul_vk)./(f_j)^(2.5);

        tmp3 = zeros(M, M);
        disp('Starting integration for equation (7.12) for ul vk')
        for i=1:M
            for j=1:M
            tmp3(i, j) = vpaintegral(internal_int3(i, j), t, 0, 2*pi);
            end   
            status = [num2str(i*100/(M)), '% completed'];
            disp(status)
        end
        Au_ul_vk = Au_ul_vk - n/4 .* tmp3;
    end
    
    Hu = sym(zeros(2*M, 2*M));
    Hu(1:M, 1:M) = Au_ul_uk;
    Hu(M+1:end, M+1:end) = Au_vl_vk;
    Hu(1:M, M+1:end) = Au_ul_vk;
    Hu(1:M , (M+1)/2) = 0;
    Hu(:, ((3*M+1)/2)) = 0;
    Hu(((3*M+1)/2), :) = 0;
    
    H = Hk + Hu;
    H(M+1:end, 1:M) = transpose(H(1:M, M+1:end));
    H = double(H);
    
    % removing rows and columns corresponding to u0 v0
    H((M+1)/2,:) = []; H(((3*M+1)/2) -1,:) = [];
    H(:, (M+1)/2) = []; H(:, ((3*M+1)/2) -1) = [];
end