function actelion_histograms
    % Software to analyse data sets using spectral parameter histograms
    % First, need to form a list of studies to examine

    actelion_filter = [2 8 15 22 26 32 38]; % IDs for the second actelion night
    control_filter = [1 2 3 4 5 6 7 8 9]; % IDs for all of the control subjects

    limits = 0:0.005:0.1; % Histogram limits. Single number is minimum spacing, vector is bin edges
    query = 'a.int_delta';

    % --------


   
    actelion = analyse_spectra('actelion',{'S3','S4'},actelion_filter,[]);
    control = analyse_spectra('actelion',{'S3','S4'},control_filter,[]);
    labels{1} = sprintf('Actelion (%i)',length(actelion.int_p));
    labels{2} = sprintf('Control (%i)',length(control.int_p));
    
    hist_overlay_query({actelion,control},labels,query,limits);
    
end

