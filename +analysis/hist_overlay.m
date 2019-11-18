function h = hist_overlay(data,f_override)
	% Make an overlay histogram of data, where data is a cell array
	% This allows different entries of data to have different amounts of data
	% f_override is a bit ugly- if it is a single number (1) then it will force
	% force a minimum value for the spacing 
	% If it is a vector, it corresponds directly to the bin edges for the histogram
	if nargin < 2 || isempty(f_override)
	    f_override = 0;
    end
    
	n = size(data,2);
	color = {'b','m','c','r','g','y','k'};

	figure
	hold on

	min_spacing = NaN;
	lower = NaN;
	upper = NaN;

	% Calculate the bin size and limits to use for each data set
	if length(f_override) == 1
	    for j = 1:n
		    [y,x] = hist(data{j},10);
		    if f_override
                min_spacing = max([min_spacing, diff(x)])
            else
                min_spacing = min([min_spacing, diff(x)]);
            end
	        lower = min([lower,x(1)]);
	        upper = max([upper,x(end)]);
	    end

	    % At most 1000 bins?
	    if (upper-lower)/1000 > min_spacing
		    min_spacing = (upper-lower)/1000;
	    end
	    bin_edges = [lower:min_spacing:upper];
	else
	    bin_edges = f_override;
    end
	
	% Make each of the histograms
	h = zeros(1,n);
	limit = [Inf -Inf];

	% First, come up with the range
	for j = 1:n
		limit(1) = min([limit(1); data{j}]);
		limit(2) = max([limit(2); data{j}]);
	end
	limit = mean(limit) + [-1 1]*1.1*diff(limit)/2;

	for j = 1:n
		this_color = color{mod(j-1,length(color))+1}; % J-1 means the first colour is yellow
	 	[y,x] = ksdensity(data{j},linspace(limit(1),limit(2),500));
	 	h(j) = plot(x,y,this_color,'LineWidth',2);
	end

	set(gca,'XLim',limit);
	yl = get(gca,'YLim');
	set(gca,'YLim',[0 yl(2)]);

