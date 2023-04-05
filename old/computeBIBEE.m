function [dG_cfa, dG_p, dG_lb] = computeBIBEE(geometry,q)
Bq = geometry.B * q;

sigma_cfa = geometry.A_cfa\Bq;
sigma_p   = geometry.A_p\Bq;
sigma_lb  = geometry.A_lb\Bq;

dG_cfa = .5 * 332.112 * q' * geometry.C * sigma_cfa;
dG_p   = .5 * 332.112 * q' * geometry.C * sigma_p;
dG_lb  = .5 * 332.112 * q' * geometry.C * sigma_lb;

