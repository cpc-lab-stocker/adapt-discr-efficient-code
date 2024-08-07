function [ y ] = sum_n_vmpdf_180( alpha, k, theta, kappa )
%sum of multiple von mises and a uniform distribution
%   y = k1*vm(kappa1) + k2*vm(kappa2) + ... + (1-sum(k))*1/180
n = length(k);
k0 = 1-sum(k);
thetapi = alpha/180*pi*2;
y = k0*1/180 * ones(size(thetapi));
for i = 1:n
    y = y + k(i)*circ_vmpdf( thetapi, theta(i)/180*pi*2, kappa(i) )/180*2*pi;
end
end

