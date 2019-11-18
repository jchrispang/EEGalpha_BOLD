function hist_overlay_query(analysis_vector,set_names,query,limits)
	% Make a overlay histogram of an arbitrary combination of parameter values
	% Take in an analysis vector
	% Evaluate the query
	% Can also specify the vector of bin edges
	if nargin < 4
	    limits = 0;
	end
	
	if isempty(set_names)
	
        % This is a bit dodgy...need to make sure these switches for set_names are correct when the data is changed
        switch length(analysis_vector)
            case 6 % control
	            set_names = {'control_W','control_R','control_S1','control_S2','control_S3','control_S4'};
	        case 5 % actelion
                set_names = {'actelion_R','actelion_S1','actelion_S2','actelion_S3','actelion_S4'};
            case 4
                set_names = {'k=0','Sum','Integrated','Int \beta=90'};
	        otherwise
	            set_names = {'220\_ec','220\_eo','1400\_ec','actelion\_R','actelion\_S1','actelion\_S2','actelion\_S3','actelion\_S4'};
        end
    end

	
	% Suppose the query was 'a.alpha./a.beta'
	% This suggests a syntax! The variable 'a' will be accessible

	data = {};
	for j = 1:length(analysis_vector) % Populate the data array
		a = analysis_vector{j};
		data{j} = eval(query);
	end
	h = hist_overlay(data,limits);
	legend(h,set_names);
	xlabel(sprintf('Parameter value (%s)',query),'interpreter','none');
	box on
end
