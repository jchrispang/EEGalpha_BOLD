function out = rho1(Q,p)
	if nargin  == 1
		p = Q;
		Q = p.phia;
	end
	
    out = -(Q.^2)./(p.sigma*p.qmax) + Q./p.sigma;
