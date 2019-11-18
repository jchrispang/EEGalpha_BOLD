function [stab,f,P,w,q2] = xyz_spectrum(self,f)
	% Formula from Roberts 2012 (Corticothalamic dynamics...)
	if nargin < 2 || isempty(f)
		f=linspace(0,100,10000);
	end
	w = f*2*pi;
	
	L2=1./((1-1i*w/self.alpha(1)).*(1-1i*w/self.beta(1))).^2; 
    Zp = self.xyz(3)*(self.alpha(1)+self.beta(1)).^2./(self.alpha(1)*self.beta(1));
    

    d = (((1-1i*w/self.gammae).^2) - self.xyz(1)).*(1+L2.*Zp) - self.xyz(2).*(1+Zp).*exp(1i*w*self.t0);
    stab=d(1)>0 && ~any(real(d(2:end))<0 & imag(d(2:end)).*imag(d(1:end-1))<0);

    q2 = ((((1-1i*w/self.gammae).^2) - self.xyz(1) -(self.xyz(2)*(1+Zp))./(1+L2.*Zp).*exp(1i*w*self.t0)))./self.re^2;

    P = 1./(abs(1+Zp.*(L2)).^2).*((1./abs(q2)).^2);
    