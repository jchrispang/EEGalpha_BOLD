function spectra = NFT_powercalc(trans_func, w, k)
% Calculate power spectra of trans_func assuming the x-axis is for w and
% the y-axis is for k

k0 = 10;

conduction = repmat(exp(-k/k0^2)', 1, length(w));

spectra.all = abs(trans_func.^2);
spectra.spatial = trapz(w, abs(trans_func.^2), 2).';
spectra.temporal = trapz(k, abs(trans_func.^2), 1)*(2*pi);
% spectra.temporal2 = trapz(k, abs(trans_func.^2).*conduction, 1)*(2*pi);