function analysis = get_spec_analysis(f_in,P_in)
	% Return a struct with all of the variables that are accessible to 
	% validate_spectrum for later analysis
    % Edited by James Pang 2017
    
	if nargin == 1 && isstruct(f_in)
		p = f_in;
		[f_in,P_in] = analytic_spectrum(p);
	end
	
	f_in = f_in(:);
	P_in = P_in(:);

    % First, reduce to a common frequency range and renormalize

    [bands,bandstr] = model.eeg_bands;

	%truncate_range = [1.36 43.36]; % Ranges in eegfit
	truncate_range = [bands(1,1) bands(end,end)]; % Upper and lower limits on EEG bands
	
% 	f = bands(1,1):0.25:bands(end,end);
    f = bands(1,1):0.1:bands(end,end);   % edited by James, increase resolution
	P = interp1(f_in,P_in,f,'linear','extrap');
	

	% f1 = find(f_in>truncate_range(1),1,'first');
	% f2 = find(f_in<truncate_range(2),1,'last');
 %    f= f_in(f1:f2);
 %    P = P_in(f1:f2);
    %P = smooth(P,3);
    analysis.int_p = trapz(f,P);
%     P = P./analysis.int_p;        % edited by James, remove normalization


	if size(P,1) > size(f,1)
		P = P';
	end
	
	delta_range = [find(f>bands(1,1),1,'first') find(f<bands(1,2),1,'last')];
	theta_range = [find(f>bands(2,1),1,'first') find(f<bands(2,2),1,'last')];
	alpha_range = [find(f>bands(3,1),1,'first') find(f<bands(3,2),1,'last')];
	sigma_range = [find(f>bands(4,1),1,'first') find(f<bands(4,2),1,'last')];
	beta_range  = [find(f>bands(5,1),1,'first') find(f<bands(5,2),1,'last')];
	gamma_range = [find(f>bands(6,1),1,'first') find(f<bands(6,2),1,'last')];

	analysis.int_delta = trapz(f(delta_range(1):delta_range(2)),P(delta_range(1):delta_range(2)));
    analysis.int_theta = trapz(f(theta_range(1):theta_range(2)),P(theta_range(1):theta_range(2)));
    analysis.int_alpha = trapz(f(alpha_range(1):alpha_range(2)),P(alpha_range(1):alpha_range(2)));
    analysis.int_sigma = trapz(f(sigma_range(1):sigma_range(2)),P(sigma_range(1):sigma_range(2)));
    analysis.int_beta  = trapz(f(beta_range(1):beta_range(2)),P(beta_range(1):beta_range(2)));
    analysis.int_gamma  = trapz(f(gamma_range(1):gamma_range(2)),P(gamma_range(1):gamma_range(2)));
    
    lin_fit = polyfit(log(f(delta_range(1):delta_range(2))),log(P(delta_range(1):delta_range(2))),1);
    analysis.lf_slope = lin_fit(1);
    
    analysis = utils.catstruct(analysis,model.get_pkf(f,P));

end
