function [A_sym G_sym, ck] = actiongradeval(n, N, w)
    syms(sym('u', [1 N]));
    u = (sym('u',[1 N]));
    syms(sym('v', [1 N]));
    v = (sym('v',[1 N]));

    k = -((N-1)/2) : ((N-1)/2);
    Ak_sym = sym('pi')*n * sum(((k+w).^2).*(u.^2 + v.^2));

    ck = [u v];

    % potential part of action

    Au_sym = 0;
    
    P = 30;
    p = 1:30;
    for j= 1:n-1                                                                                                                                
        syms t real
        a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
        b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);
        f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
        tmp = (f_j)^(-0.5);

        % use trapezoidal rule instead of integration for P = 30
        trapz_f_j = subs(tmp, t, 2*pi.*p./P);
        integral = sum(trapz_f_j);

        % summation over j
        Au_sym = Au_sym + n*pi/P * integral;
    end 

    % combining kinetic and potential part together 
    A_sym = Ak_sym + Au_sym;
    
    % grad 
    % grad kinetic part
    Ak_uk = 2*sym('pi')*n.*((k+w).^2).*u; 
    Ak_vk = 2*sym('pi')*n.*((k+w).^2).*v;
    grad_Ak = [Ak_uk Ak_vk];

    %grad potential part
    grad_Ap = sym(zeros(2*N, 1));

    for j=1:n-1 % summation over j
        k = -((N-1)/2) : ((N-1)/2);
        syms t real
        a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
        b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);
        f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
        power_f_j = f_j^(1.5);

        %grad_f_j
        l = -((N-1)/2) : ((N-1)/2);
        grad_f_j = sym(zeros(2*N, 1));
        
        for k=-((N-1)/2) : ((N-1)/2)
            syms t real
            a_lj = (1-cos(2*sym('pi')*j.*l./n)).*cos(l.*t) + sin(2*sym('pi').*j.*l./n).*sin(l.*t);
            b_lj = (-1+cos(2*sym('pi')*j.*l./n)).*sin(l.*t) + sin(2*sym('pi').*j.*l./n).*cos(l.*t);
            a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
            b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);

            tmp = 2 * sum((a_lj.*a_kj + b_lj.*b_kj).*u + (b_lj.*a_kj - a_lj.*b_kj).*v);
            grad_f_j(k+(N+1)/2) = tmp;
            
            tmp2 = 2 * sum((a_lj.*b_kj - b_lj.*a_kj).*u + (b_lj.*b_kj + a_lj.*a_kj).*v);
            grad_f_j(k+(3*N+1)/2) = tmp2;
        end

        expression_inside = grad_f_j./power_f_j;

        disp('Evaluating integral for gradient')
        trapz_integral = subs(expression_inside, t, 2*pi.*p./P);
        integral_part = sum(trapz_integral, 2);
        
        %for i=1:(2*N)
            %integral_part(i) = vpaintegral(expression_inside(i), t, 0, 2*pi);
            %status = [num2str(i),'/',num2str(2*N), ' evaluated'];
            %disp(status)
        %end
        disp('Done')
        final_expr_j = -((n*pi)/(2*P)).*integral_part;
        grad_Ap = grad_Ap + final_expr_j; %summation over j via loop
    end
    G_sym = grad_Ak + transpose(grad_Ap);  
end