clear all; clc; % clean the memory and the screen

% define the exact solution for verification purposes
kappa = @(x) 1 + x.^2;
exact = @(x) x.^4;
exact_x = @(x) 4*x.^3;
%exact = @(x) x.^3;
%exact_x = @(x) 3*x.^2;

% define the loading and boundary info
%f = @(x) -6*x;        % u,xx + f = 0
f = @(x) -12*x.^2 - 20*x.^4;
g = exact(1);           %     u(1) = g
h = -exact_x(0);        %  -u,x(0) = h

% quadrature rule
n_int = 3;
[xi, weight] = Gauss(n_int, -1, 1);

% interpolation degree
pp   = 1;                     % degree of shape function
n_el = 5;                     % number of elements
n_np = pp * n_el + 1;         % number of nodal points
n_en = pp + 1;                % number of element nodes
n_eq = n_np - 1;              % number of equations
x_coor = 0 : 1/(n_el*pp) : 1; % nodal coordinates

IEN = zeros(n_en, n_el);

for aa = 1 : n_en
    for ee = 1 : n_el
        IEN(aa,ee) = (ee-1)*pp + aa;
    end
end

ID = 1 : n_np;
ID(end) = 0;

% allocate zero matrix and zero vector
K = spalloc(n_eq, n_eq, min(n_eq^2, (2*pp + 1)*n_eq));
F = zeros(n_eq, 1);

for ee = 1 : n_el
    k_ele = zeros(n_en, n_en);
    f_ele = zeros(n_en, 1);
    x_ele(:) = x_coor(IEN(:,ee));

    for ll = 1 : n_int
        dx_dxi = 0.0;
        x_qua  = 0.0;
        for aa = 1 : n_en
            x_qua  = x_qua  + x_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : n_en
            f_ele(aa) = f_ele(aa) + weight(ll)*kappa(x_qua) * PolyShape(pp, aa, xi(ll), 0) * ...
                f(x_qua) * dx_dxi;

            for bb = 1 : n_en
                k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * ...
                    PolyShape(pp, aa, xi(ll), 1) * PolyShape(pp, bb, xi(ll), 1) *...
                    dxi_dx;
            end
        end
    end  % end of quadrature

    for aa = 1 : n_en
        PP = ID(IEN(aa, ee));
        if PP > 0
            F(PP) = F(PP) + f_ele(aa);

            for bb = 1 : n_en
                QQ = ID(IEN(bb, ee));
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_ele(aa,bb);
                else
                    F(PP) = F(PP) - k_ele(aa,bb)*g;
                end
            end
        end
    end
end 

F(1) = F(1) + h;

uh_temp = K \ F; % solve K uh = F linear system

uh = [uh_temp; g];

% calculate relative L2, H1 seminorm, and full H1 errors
n_int = pp+4;
[xi, weight] = Gauss(n_int, -1, 1);

errL2sq = 0.0; normL2sq = 0.0;
errH1semiSq = 0.0; normH1semiSq = 0.0;
for ee = 1 : n_el
    x_ele(:) = x_coor(IEN(:,ee));
    u_ele(:) = uh(IEN(:,ee));

    for ll = 1 : n_int
        x = 0.0; dx_dxi = 0.0;
        uh_q = 0.0; duh_dxi = 0.0;
        for aa = 1 : n_en
            N = PolyShape(pp, aa, xi(ll), 0);
            dN_dxi = PolyShape(pp, aa, xi(ll), 1);
            x = x + x_ele(aa) * N;
            dx_dxi = dx_dxi + x_ele(aa) * dN_dxi;
            uh_q = uh_q + u_ele(aa) * N;
            duh_dxi = duh_dxi + u_ele(aa) * dN_dxi;
        end
        duh_q = duh_dxi / dx_dxi;
        dV = weight(ll) * dx_dxi;
        u_exact = exact(x);
        du_exact = exact_x(x);

        errL2sq = errL2sq + (uh_q - u_exact)^2 * dV;
        normL2sq = normL2sq + u_exact^2 * dV;
        errH1semiSq = errH1semiSq + (duh_q - du_exact)^2 * dV;
        normH1semiSq = normH1semiSq + du_exact^2 * dV;
    end
end

% Relative errors require nonzero exact-solution norms.
assert(normL2sq > 0, 'Relative L2 error requires a nonzero exact-solution norm.');
assert(normH1semiSq > 0, ...
    'Relative H1 seminorm error requires a nonzero exact-derivative norm.');

eL2 = sqrt(errL2sq / normL2sq);
eH1semi = sqrt(errH1semiSq / normH1semiSq);
eH1 = sqrt((errL2sq + errH1semiSq) / (normL2sq + normH1semiSq));

fprintf('Relative L2 error: %.6e\n', eL2);
fprintf('Relative H1 seminorm error: %.6e\n', eH1semi);
fprintf('Relative H1 error: %.6e\n', eH1);

% EOF
