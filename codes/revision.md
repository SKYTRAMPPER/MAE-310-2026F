1. 你至少需要一个pdf文档结合图片和文字来回答后三个小问
2. driver2a.m中
for aa = 1 : n_en
            f_ele(aa) = f_ele(aa) + weight(ll)*kappa(x_qua) * PolyShape(pp, aa, xi(ll), 0) * ...
                f(x_qua) * dx_dxi;

            for bb = 1 : n_en
                k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * ...
                    PolyShape(pp, aa, xi(ll), 1) * PolyShape(pp, bb, xi(ll), 1) *...
                    dxi_dx;
            end
        end
kappa 出现在了错误的位置