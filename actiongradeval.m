function objfun = actiongradeval(n, N, w)
        % kinetic part of action
    syms(sym('u', [1 N]));
    u = (sym('u',[1 N]));
    syms(sym('v', [1 N]));
    v = (sym('v',[1 N]));

    k = -((N-1)/2) : ((N-1)/2);
    Ak_sym = sym('pi')*n * sum(k.*(u.^2 + v.^2));

    x = transpose([u,v]);

    % verifying that this is indeed agrees with the explicit representation in
    % the book
    gradf = jacobian(Ak_sym,x).';

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
    %A = matlabFunction(A_sym, 'vars', {[u v] t}); %converting this into an anonymous function to be used for fminunc

    % for gradient

    % grad 
    % grad kinetic part
    Ak_uk = 2*sym('pi')*n.*((k+w).^2).*u; 
    Ak_vk = 2*sym('pi')*n.*((k+w).^2).*v;
    grad_Ak = [Ak_uk Ak_vk];

    %grad potential part
    k = -((N-1)/2) : ((N-1)/2);
    grad_Ap = sym(zeros(2*N, 1));

    for j=1:n-1 % summation over j

        syms t real
        a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
        b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);
        f_j = (sum(a_kj .* u + b_kj .* v))^2 + (sum(a_kj .* v - b_kj .* u))^2;
        power_f_j = f_j^(1.5);

        %grad_f_j
        l = -((N-1)/2) : ((N-1)/2);
        grad_f_j = sym(zeros(2*N, 1));
        disp('First')
        for k=-((N-1)/2) : ((N-1)/2)
            syms t real
            a_lj = (1-cos(2*sym('pi')*j.*l./n)).*cos(l.*t) + sin(2*sym('pi').*j.*l./n).*sin(l.*t);
            b_lj = (-1+cos(2*sym('pi')*j.*l./n)).*sin(l.*t) + sin(2*sym('pi').*j.*l./n).*cos(l.*t);
            a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
            b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);

            tmp = 2 * sum((a_lj.*a_kj + b_lj.*b_kj).*u + (b_lj.*a_kj - a_lj.*b_kj).*v);
            grad_f_j(k+(N+1)/2) = tmp;
            disp('Yay')
        end

        for k=-((N-1)/2) : ((N-1)/2)
            syms t real
            a_lj = (1-cos(2*sym('pi')*j.*l./n)).*cos(l.*t) + sin(2*sym('pi').*j.*l./n).*sin(l.*t);
            b_lj = (-1+cos(2*sym('pi')*j.*l./n)).*sin(l.*t) + sin(2*sym('pi').*j.*l./n).*cos(l.*t);
            a_kj = (1-cos(2*sym('pi')*j.*k./n)).*cos(k.*t) + sin(2*sym('pi').*j.*k./n).*sin(k.*t);
            b_kj = (-1+cos(2*sym('pi')*j.*k./n)).*sin(k.*t) + sin(2*sym('pi').*j.*k./n).*cos(k.*t);

            tmp = 2 * sum((a_lj.*b_kj - b_lj.*a_kj).*u + (b_lj.*b_kj - a_lj.*a_kj).*v);
            grad_f_j(k+(3*N+1)/2) = tmp;
            disp('Yay2')
        end

        expression_inside = grad_f_j./power_f_j; % transposed grad fj

        integral_part = int(expression_inside, t, 0, 2*pi); 
        final_expr_j = -(n/4).*integral_part;
        grad_Ap = grad_Ap + final_expr_j; %summation over j via loop
        disp('Second')
    end

    G_sym = transpose(grad_Ak) + grad_Ap;

    objfun = matlabFunction(A_sym, G_sym,'vars',{[u v] t},'Outputs',{'f','gradf'});
    
end