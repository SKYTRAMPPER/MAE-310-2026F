clear all; clc; close all;

kappa   = @(x) 1 + x.^2;
exact   = @(x) x.^4;
exact_x = @(x) 4*x.^3;
f       = @(x) -12*x.^2 - 20*x.^4;
g       = exact(1);
h       = -kappa(0)*exact_x(0);

n_int_err = 20;
[xi_e, w_e] = Gauss(n_int_err, -1, 1);

pp   = 3;   n_el = 8;   n_en = pp + 1;
n_np = pp*n_el + 1;   n_eq = n_np - 1;
x_coor = linspace(0, 1, n_np);
IEN = zeros(n_en, n_el);
for aa = 1:n_en
  for ee = 1:n_el
    IEN(aa,ee) = (ee-1)*pp + aa;
  end
end
ID = 1:n_np;   ID(end) = 0;

nqs = 1:6;
fprintf('\n(d)  p=3, n_el=8, 误差积分固定 20 点\n');
fprintf(' n_q  rank(k_ele)   rcond(K)      e_L2        e_H1semi     解\n');

for m = 1:numel(nqs)
  nq = nqs(m);
  [xi, weight] = Gauss(nq, -1, 1);

  K = spalloc(n_eq, n_eq, min(n_eq^2, (2*pp+1)*n_eq));
  F = zeros(n_eq, 1);
  k_first = zeros(n_en, n_en);

  for ee = 1:n_el
    k_ele = zeros(n_en, n_en);
    f_ele = zeros(n_en, 1);
    x_ele = x_coor(IEN(:,ee));

    for ll = 1:nq
      x_qua = 0.0;  dx_dxi = 0.0;
      for aa = 1:n_en
        x_qua  = x_qua  + x_ele(aa)*PolyShape(pp,aa,xi(ll),0);
        dx_dxi = dx_dxi + x_ele(aa)*PolyShape(pp,aa,xi(ll),1);
      end
      dxi_dx = 1.0/dx_dxi;
      for aa = 1:n_en
        f_ele(aa) = f_ele(aa) + weight(ll)*PolyShape(pp,aa,xi(ll),0)*...
            f(x_qua)*dx_dxi;
        for bb = 1:n_en
          k_ele(aa,bb) = k_ele(aa,bb) + weight(ll)*kappa(x_qua)*...
              PolyShape(pp,aa,xi(ll),1)*PolyShape(pp,bb,xi(ll),1)*dxi_dx;
        end
      end
    end

    for aa = 1:n_en
      PP = ID(IEN(aa,ee));
      if PP > 0
        F(PP) = F(PP) + f_ele(aa);
        for bb = 1:n_en
          QQ = ID(IEN(bb,ee));
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

  rk = rank(full(k_first));
  rc = rcond(full(K));
  try
    uh_temp = K \ F;
  catch
    uh_temp = nan(n_eq,1);
  end
  uh = [uh_temp; g];
  ok = all(isfinite(uh));
  eL2 = NaN;  eH1 = NaN;
  if ok
    a=0; b=0; c=0; d=0;
    for ee = 1:n_el
      x_ele = x_coor(IEN(:,ee));
      u_ele = uh(IEN(:,ee));
      for ll = 1:n_int_err
        x=0; dx_dxi=0; uh_q=0; duh_dxi=0;
        for aa = 1:n_en
          N  = PolyShape(pp,aa,xi_e(ll),0);
          dN = PolyShape(pp,aa,xi_e(ll),1);
          x       = x       + x_ele(aa)*N;
          dx_dxi  = dx_dxi  + x_ele(aa)*dN;
          uh_q    = uh_q    + u_ele(aa)*N;
          duh_dxi = duh_dxi + u_ele(aa)*dN;
        end
        duh_q = duh_dxi/dx_dxi;
        dV = w_e(ll)*dx_dxi;
        a = a + (uh_q - exact(x))^2*dV;
        b = b + exact(x)^2*dV;
        c = c + (duh_q - exact_x(x))^2*dV;
        d = d + exact_x(x)^2*dV;
      end
    end
    eL2 = sqrt(a/b);   eH1 = sqrt(c/d);
    fprintf('%4d      %d       %10.2e   %10.4e  %10.4e   OK\n', nq, rk, rc, eL2, eH1);
  else
    fprintf('%4d      %d       %10.2e        --          --      Inf/NaN (K singular)\n', ...
            nq, rk, rc);
  end
end

% EOF
