function T = NFT_transfunc(p, w, k)
% Calculate NFT transfer function.

[ww, kk] = meshgrid(w, k);

p.Gesn = 6.16*0.4323; % default value
p.Gee = p.gabcd(1);
p.Gei = p.gabcd(2);
p.Gese = p.gabcd(3);
p.Gesre = p.gabcd(4);
p.Gsrs = p.gabcd(5);

L = (1 - 1i*ww/p.alpha(1)).^(-1).*(1 - 1i*ww/p.beta(1)).^(-1);
A = (p.Gesn).*(L.^2).*exp(1i*ww*p.t0/2)./  ...
    ((1 - p.Gsrs*L.^2).*(1 - p.Gei*L));
q2 = (1 - 1i*ww/p.gammae).^2 - ...
    (1./(1 - p.Gei*L)).*(p.Gee*L + ...
           (p.Gese*L.^2 + p.Gesre*L.^3).*exp(1i*ww*p.t0)./ ...
           (1 - p.Gsrs*L.^2));

T = A./(kk.^2 + q2);

