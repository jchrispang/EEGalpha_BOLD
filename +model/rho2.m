function out = rho2(Q,p)
    out = 2*Q.^3./p.sigma^2./p.qmax^2 - 3*Q.^2./p.sigma^2/p.qmax + Q/p.sigma^2;
