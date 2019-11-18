% This is a USEFUL snippet which generates statistical analyses for
% each of the experimental data sets (e.g. Actelion)
sets = {'220_ec','220_eo','1400_ec','actelion_R','actelion_S1','actelion_S2','actelion_S3','actelion_S4'};
sets = {'actelion_R','actelion_S1','actelion_S2','actelion_S3','actelion_S4'};
sets = {'control_W','control_R','control_S1','control_S2','control_S3','control_S4'};
sets = {'br_EC','br_EO','control_S2','control_S3'};
sets = {'br_EC','br_EO','control_R','control_S1','control_S2',{'control_S3','control_S4'}};

for j = 1:length(sets)
	analyse{j} = analyse_spectra(sets{j});
	%print_analysis(analyse{j});
end
%make_hist_graphs(analyse)

%No scaling parameter for the limits. Scaling is meaningless for the analytic spectrum
attribs = {'deltamax','thetamax','alphamax','sigmamax','betamax','int_delta','int_theta','int_alpha','int_sigma','int_beta','delta_str','theta_str','alpha_str','sigma_str','beta_str','int_gamma','lf_slope','delta_enhance'};
for j = 1:length(analyse)
	sets{j}
	for k = 1:length(attribs)
		vals = prctile(analyse{j}.(attribs{k}),[10 90]);
		if strfind(attribs{k},'str')
			fprintf('( ~isfinite(a.%s) | (a.%s >= %.2f & a.%s <= %.2f));...\n',attribs{k},attribs{k},vals(1),attribs{k},vals(2));
		elseif strfind(attribs{k},'max')
			fprintf('(1 | a.%s >= %.2f & a.%s <= %.2f);...\n',attribs{k},vals(1),attribs{k},vals(2));
		else
			fprintf('(a.%s >= %.2f & a.%s <= %.2f);...\n',attribs{k},vals(1),attribs{k},vals(2));
		end
	end
end
