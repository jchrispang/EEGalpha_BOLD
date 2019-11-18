function fig_handles = make_hist_graphs(analysis_vector,set_names,attribs,titles)
	% make_hist_graphs(analysis_vector,set_names,attribs,titles,prefix)
	% analysis_vector is a cell array of analysis objects
	% set_names is a cell array the same size as analysis_vector with the set names to use for the legend
	% attribs is a cell array of all the Matlab variables plots are required for e.g. int_delta
	% titles is a cell array of the names of those variables that appear on the x-axis label
	% Internal definitions of attribs and titles in Rev 852
	if nargin < 4 || isempty(titles)
		titles = {};
	end
	
	% All other arguments are mandatory
	
	for k = 1:length(attribs)
		data = {};
		for j = 1:length(analysis_vector) % Populate the data array
			data{j} = analysis_vector{j}.(attribs{k});
		end

		if any(strfind(attribs{k},'max'))
		        spacing = [0:0.05:1]; % To force analytic data to display on comparable axes
		        spacing = 0; % For experimental data
		elseif any(strfind('f_',attribs{k}));
		    spacing = 1;
		else
		    spacing = 0;
		end
		    
		h = hist_overlay(data,spacing);

		legend(h,set_names,'Interpreter','none');
		%title(sprintf('%s %s',titles{k},strrep(prefix,'_','\_')));
		%xlabel(sprintf('Parameter value (%s) ',attribs{k}));
		xlabel(strrep(titles{k},'$',''));
		set(gca,'FontSize',20)
		set(gcf,'Renderer','Painters');
		box on
		
		% Code to check ANOVAs
		%{
		anova_input = pad_nan(data([1,2])); % Select which columns to compare
	    %anova_input = pad_nan(data()); % All columns
		fprintf('\nStatistic: %s, ANOVA p=%f\n',titles{k},anova1(anova_input,[],'on'))
        set(figure(1),'Position',[ 9   379   859   670]);
        set(figure(2),'Position',[ 910   666   560   420]);
        set(figure(3),'Position',[ 909   140   560   420]);
		pause
		%}
		%anova_input = pad_nan(data([1,2])); % Select which columns to compare- here, only first two
		%anova_input = pad_nan(data());
		%p_value = anova1(anova_input,[],'off');
		%text(0.9,0.5,sprintf('ANOVA P=%g',p_value),'Units','normalized','HorizontalAlignment','right');
		
		pfig(sprintf('kde_%i',k));
		%saveas(gcf,sprintf('/home/romesha/Desktop/histarea_%s_%s_%i_2.fig',prefix,date,k));
		close all
	end	
end
	
function m = pad_nan(data)
	% First, get max number of entries
	n = 0;
	for j = 1:length(data)
	    n = max(n,length(data{j}));
	    %n = min(n,length(data{j}));
    end
    % Write it to a matrix
    m = NaN(n,length(data));
    
    for j = 1:length(data)
        m(:,j) = [data{j};NaN(n-length(data{j}),1)];
        %a = data{j};
        %m(:,j) = a(1:n);
    end
end
    

