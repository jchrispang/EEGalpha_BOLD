function plot_containers = cloud_corr(data,lim,plot_type,binres,colour,titles)
	% Plot a correlation matrix of cloudplots (column-wise)
	%
	% cloud_corr(data,uselog,plot_type,binres,colour)
	%
	% uselog = 1 for a log scale of density
	% plot_type: this is 1 for nuab, 2 for gab. Tries to autodetect based on mean value of first column if not specified
	% binres: override the default resolution (default 250x250px per plot)
	% colour: override the plot with some other colour axis (e.g. phi)
	%
	% This code *should* be able to overlay plots directly
	% However, nasty things will probably happen if the binres changes between plots
	% or if the axis limits changes, or if a non-binary colour is used i.e. base plots have phi

	uselog = 0; % Don't use log by default, normally 'canvas = canvas>0' is set in cloud_plots.m

	if isstruct(data) % If a plot_containers object is given instead of a matrix
		mdim = sqrt(length(data));
		if nargin < 3
		    plot_type = log(mean(data(1).axisLimits)) > 0; % Autodetect if not specified
		end
	else
		mdim = size(data,2)-1;
		if nargin < 3 || isempty(plot_type)
		    plot_type = log(mean(data(:,1))) > 0; % Autodetect if not specified
		end
	end
	
	if nargin < 5 || isempty(colour)
		colour = [];
	end
	
	if nargin < 4 || isempty(binres)
		binres = 100; % assume square
	end

	if nargin == 6 % If titles are specified, automatic limit
	    lim = []; 
	else % Otherwise, automatically determine titles and limits
		if mdim == 3 % phi
				titles = {'\Phi_e','\Phi_r','\Phi_s'};
				if nargin < 2 || isempty(lim)
					lim = limits(3);
				end
	    elseif mdim == 5 % gabcd
	    		titles = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}'};
				if nargin < 2 || isempty(lim)
					lim = [-20 -20 -20  -20 -20;...
					       20   20  20   20  20 ];
				end
		elseif plot_type == 0 % nuab 
				titles = {'\nu_{ee}','\nu_{ei}','\nu_{es}','\nu_{se}','\nu_{sr}','\nu_{sn}','\nu_{re}','\nu_{rs}'};
				if nargin < 2 || isempty(lim)
					lim = limits(2);
				end
		elseif plot_type == 1 % gab 
				titles = {'G_{ee}','G_{ei}','G_{es}','G_{se}','G_{sr}','G_{sn}','G_{re}','G_{rs}'};
				if nargin < 2 || isempty(lim)
					lim = limits(1);
				end
		elseif nargin < 6 || isempty(titles)
			error('Could not determine plot type')
		end
	end

	if isstruct(data) % If we are simply showing an existing plot
		plot_containers = show_plots(data,titles);
		return
	end
	

	% Preallocate the container
	plot_containers(mdim^2).canvas = [];
	plot_containers(end).axisLimits = [];
	plot_containers(end).colormap = [];
	
	parfor j = 1:mdim^2 % For each bin index
		[row,col] = get_dims(j,mdim);
		
		if col < row
			if ~isempty(lim)
				axis_limits = [lim(1,col),lim(2,col),lim(1,row),lim(2,row)];
			else
				axis_limits = [];
			end		
			plot_containers(j) = cloud_plots(data(:,col),data(:,row),axis_limits,uselog,[binres binres],colour,1,1);
		end
	end
	
	show_plots(plot_containers,titles);
	
end

function [row,col] = get_dims(j,mdim)
		if rem(j,mdim) == 0
			col = mdim;
			row = (j/mdim);
		else
			col = rem(j,mdim);
			row = floor(j/mdim)+1;
		end
		row = row+1;
end

function plot_containers = show_plots(plot_containers,titles)
	nplots = length(findall(gcf,'type','axes')); % NUMBER OF SUBPLOTS
	binres = size(plot_containers(1).canvas,1);
	mdim = sqrt(length(plot_containers));
	
	if nplots > 1 % There are existing plots- Do some tricky stuff!
		if length(plot_containers) ~= nplots; error('Wrong number of plots to overlay!'); end

		% First, we need to regenerate plot_containers based on the existing plot
		existing_containers(mdim^2).canvas = [];
		existing_containers(end).axisLimits = [];
		
		% Is this the first plot?
		subplot(mdim,mdim,1);
		c = get(gca,'Children');
		canvas = get(c(end),'CData');
		first_overlay = max(canvas(:))==1; % Is the maximum value 1?
		% This determines whether we need to invert the canvas or not (imcomplement)
		% The reason for this is because the first plot is a binary imcomplement from cloud_plots
		% After an overlay is done, there is colour data so inverting is no longer necessary
		
		for j = 1:mdim^2
			[row,col] = get_dims(j,mdim);
			subplot(mdim,mdim,j)
			if col < row
				c = get(gca,'Children');
				canvas =  rot90(fliplr(get(c(end),'CData')));
				if first_overlay
				    existing_containers(j).canvas = imcomplement(canvas);
				else
                    existing_containers(j).canvas = canvas;
                end
				existing_containers(j).axisLimits = [get(gca,'XLim'),get(gca,'YLim')];
			end
		end
		
		% Assume that the existing data is in binary form i.e. 1,2,4,8 etc. max
		n_existing = get(gca,'UserData'); % Number of existing plots
		multiplier = 2^(n_existing); % Next power of 2
		
		% So we know that colours will form a continuum of powers of 2. Let's just go ahead and define some

		for j = 1:mdim^2
			% Calculate the new plot
			plot_containers(j).canvas = imcomplement(plot_containers(j).canvas)*multiplier + existing_containers(j).canvas;
			plot_containers(j).colormap = jet; % Switdch to colourmap
		end
		n_data = n_existing + 1;

	else
		n_data = 1; % Number of data series on the plot
	end
		
	set(gcf,'Position',[1 1 binres*mdim binres*mdim]);
	clf
	
	color_range = [min(min(min([plot_containers.canvas]))) max(max(max([plot_containers.canvas]))) ];
	color_range = color_range + [-diff(color_range) diff(color_range)]*0.1;
	
	for j = 1:mdim^2 % For each bin index
		[row,col] = get_dims(j,mdim);
		subplot(mdim,mdim,j)
		if col >= row
			set(gca,'Visible','off')
		else
			cloud_plots(plot_containers(j));
			set(gca,'CLim',color_range);
			add_title(row,col,mdim,titles);
			set(gca,'UserData',n_data);
		end

	end
	
	set(gcf,'Position',[1 1 binres*mdim binres*mdim]) % 250 pixels per plot
end

function add_title(row,col,mdim,titles)
	if row == mdim+1 % If this is the last row, label goes on the bottom
		text('units','normalized','pos',[0.5,-0.5],'string',titles{col},'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
	end
	
	if col == 1
		text('units','normalized','pos',[-0.4 0.5],'string',titles{row},'FontSize',16,'HorizontalAlignment','right','VerticalAlignment','middle')
	end
	
end
