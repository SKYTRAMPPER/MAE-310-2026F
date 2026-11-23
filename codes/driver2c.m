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
n_int_err = 20;
[xi_e, w_e] = Gauss(n_int_err, -1, 1);
n_els = [2 4 6 8 10 12 14 16];
hh = 1 ./ n_els;
pps   = [1 2 3];
errL2 = zeros(numel(pps), numel(n_els));
errH1 = zeros(numel(pps), numel(n_els));
slopeL2 = zeros(numel(pps),1);  slopeH1 = zeros(numel(pps),1);
mkr = {'o','s','^'};
figure; hold on; grid on;
mk  = {'-o', '-s', '-^'};
col = {[0.11 0.36 0.65], [0.13 0.60 0.47], [0.80 0.30 0.20]};

% interpolation degree
%n_eq = n_np - 1;              % number of equations
%x_coor = 0 : 1/(n_el*pp) : 1; % nodal coordinates

for ip = 1:numel(pps)
  pp   = pps(ip);
  n_en = pp + 1;
  n_int_asm = pp + 4;
  [xi, weight] = Gauss(n_int_asm, -1, 1);

  for k = 1:numel(n_els)
    n_el = n_els(k);
    n_np   = pp * n_el + 1;
    n_eq   = n_np - 1;
    x_coor = linspace(0, 1, n_np);
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
    x_ele = x_coor(IEN(:,ee));

for ll = 1 : n_int_asm
        dx_dxi = 0.0;
        x_qua  = 0.0;
        for aa = 1 : n_en
            x_qua  = x_qua  + x_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : n_en
            f_ele(aa) = f_ele(aa) + weight(ll) * PolyShape(pp, aa, xi(ll), 0) * ...
                f(x_qua) * dx_dxi;

            for bb = 1 : n_en
                k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * kappa(x_qua) * ...
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
errL2sq = 0.0; normL2sq = 0.0;
errH1semiSq = 0.0; normH1semiSq = 0.0;
for ee = 1 : n_el
    x_ele = x_coor(IEN(:,ee));
    u_ele = uh(IEN(:,ee));

    for ll = 1 : n_int_err
        x = 0.0; dx_dxi = 0.0;
        uh_q = 0.0; duh_dxi = 0.0;
        for aa = 1 : n_en
            N = PolyShape(pp, aa, xi_e(ll), 0);
            dN_dxi = PolyShape(pp, aa, xi_e(ll), 1);
            x = x + x_ele(aa) * N;
            dx_dxi = dx_dxi + x_ele(aa) * dN_dxi;
            uh_q = uh_q + u_ele(aa) * N;
            duh_dxi = duh_dxi + u_ele(aa) * dN_dxi;
        end
        duh_q = duh_dxi / dx_dxi;
        dV = w_e(ll) * dx_dxi;
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
assert(normH1semiSq > 0, 'Relative H1 seminorm error requires a nonzero exact-derivative norm.');
eL2     = sqrt(errL2sq / normL2sq);
eH1semi = sqrt(errH1semiSq / normH1semiSq);
errL2(ip,k) = eL2;
errH1(ip,k) = eH1semi;
  end
  loglog(hh, errL2(ip,:), ['-' mkr{ip}], 'Color', col{ip}, 'LineWidth', 1.4, ...
         'MarkerSize', 6, 'DisplayName', sprintf('p=%d  e_{L2}', pp));
  loglog(hh, errH1(ip,:), '--', 'Color', col{ip}, 'LineWidth', 1.4, ...
         'Marker', mkr{ip}, 'MarkerSize', 5, 'DisplayName', sprintf('p=%d  e_{H1}', pp));
  cL2 = polyfit(log(hh), log(errL2(ip,:)), 1);
  cH1 = polyfit(log(hh), log(errH1(ip,:)), 1);
  slopeL2(ip) = cL2(1);   slopeH1(ip) = cH1(1);
  fprintf('p=%d  slope e_L2 = %.3f (theory %d),  slope e_H1 = %.3f (theory %d)\n', ...
          pp, cL2(1), pp+1, cH1(1), pp);
end   % end of for ip（p = 1, 2, 3）

xlabel('mesh size h');  ylabel('relative error');
title('HW3 (b)(c): convergence for p = 1, 2, 3');
legend('Location','southoutside','NumColumns',2,'FontSize',8);
grid on;

fprintf('\n=== relative errors (each p, 8 meshes) ===\n');
for ip = 1:numel(pps)
  fprintf('\n p = %d\n', pps(ip));
  fprintf('  n_el      h         e_L2        e_H1semi\n');
  fprintf('%5d  %8.5f  %11.4e  %11.4e\n', ...
          [n_els; hh; errL2(ip,:); errH1(ip,:)]);
end

fprintf('\n=== convergence rates ===\n');
fprintf('   p   slope L2   theory(p+1)   slope H1   theory(p)\n');
for ip = 1:numel(pps)
  fprintf('%4d   %8.3f   %8d      %8.3f   %8d\n', ...
          pps(ip), slopeL2(ip), pps(ip)+1, slopeH1(ip), pps(ip));
end


% EOF
