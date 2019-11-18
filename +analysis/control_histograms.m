function snippet8()
    % This snippet constructs histograms for the control data
    % SVN rev 223 had the old functionality (including fitting refits)

    fit_sets = {'control_S3_fitted','control_S3_fitted_k0_20','control_S3_fitted_k0_15','control_S3_fitted_k0_10','control_S3_fitted_k0_5','control_S3_fitted_k0_1'};

    % Get analysis
    control_fits = {};

    for j = 1:length(fit_sets)
        fprintf(1,'Analysing %s (as fitted)\n',fit_sets{j});
        control_fits{j} = get_analysis_from_fit(fit_sets{j});
    end
    
    % Need to use the same frequency limits for the RAW spectrum, or else the powers don't match
    baseline = get_exp_analysis_from_fit('control_S3_fitted');
    
    % Make plots
    make_hist_graphs([{baseline} control_fits],[{'control_S3'},fit_sets]);
end

function out = convert_analysis(in)
    attribs = {'deltamax','thetamax','alphamax','sigmamax','betamax','f_delta','f_theta','f_alpha','f_sigma','f_beta','int_delta','int_theta','int_alpha','int_sigma','int_beta','int_gamma','int_p','lf_slope','hf_slope','delta_enhance'};

    for x = 1:length(attribs)
        out.(attribs{x}) = readstruct(in,attribs{x})';
    end
end

function analyse = get_analysis_from_fit(set_str)
     % Take in a string specifying which control fits to load
     % Examine the spectrum from each of them
     % Return analysis based on fits.fit
    fit_file = load(sprintf('./data/%s',set_str));
    fits = fit_file.fits;
    
    for j = 1:length(fits)
        if ~isempty(fits(j).Freq)
            f = fits(j);
            P = f.Fit/mex_trapz(f.Freq,f.Fit);
            analyse(j) = get_spec_analysis(f.Freq,P);
        end
    end
    analyse = convert_analysis(analyse);
end

function analyse = get_exp_analysis_from_fit(set_str)
     % Take in a string specifying which control fits to load
     % Examine the spectrum from each of them
     % Return analysis based on fits.fit
    fit_file = load(sprintf('./data/%s',set_str));
    fits = fit_file.fits;
    
    for j = 1:length(fits)
        if ~isempty(fits(j).Freq)
            f = fits(j);
            P = f.Exp/mex_trapz(f.Freq,f.Exp);
            analyse(j) = get_spec_analysis(f.Freq,P);
        end
    end
    analyse = convert_analysis(analyse);
end

