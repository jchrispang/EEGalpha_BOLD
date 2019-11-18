function T = get_transfunc(data, w, k)
% Calculate NFT transfer function from Brain Resource data

[ww, kk] = meshgrid(w, k);

Gesn = 6.16*0.4323; % default value
gammae = 116;       % default value

L = (1 - 1i*ww/data.Alpha).^(-1).*(1 - 1i*ww/data.Beta).^(-1);
A = (Gesn).*(L.^2).*exp(1i*ww*data.t0/2)./  ...
    ((1 - data.Gsrs*L.^2).*(1 - data.Gei*L));
q2 = (1 - 1i*ww/gammae).^2 - ...
    (1./(1 - data.Gei*L)).*(data.Gee*L + ...
           (data.Gese*L.^2 + data.Gesre*L.^3).*exp(1i*ww*data.t0)./ ...
           (1 - data.Gsrs*L.^2));

T = A./(kk.^2 + q2);