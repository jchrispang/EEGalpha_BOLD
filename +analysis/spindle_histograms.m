function spec_analysis = spindle_histograms()
	% This function produces an analysis struct suitable for histograms
	% corresponding to the sleep spindle spectra

	for j = 1:9
		[s_f,fftx_out,P_out] = spindle_readevents('control_apnea',j);
		s_P = mean(P_out,2);
		analysis(j) = get_spec_analysis(s_f,s_P);
	end
	spec_analysis = analysis;
	attribs = [fields(spec_analysis(1))];
	for x = 1:length(attribs)
		a.(attribs{x}) = readstruct(spec_analysis,attribs{x})';
	end
	spec_analysis = a;

