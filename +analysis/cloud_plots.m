function plotdata = cloud_plots(X,Y,axisLimits,useLogScale,bins,colour,suppress,bw)
	% Display a cloud plot from X and Y data
	%
	% cloudplots(X,Y,axisLimits,useLogScale,bins,colour,suppress)
	%
	% MODES- Show existing, make cloudplot, make cloudplot with colour data
	% Colour data is usually some vector of phi values (same size as Y)
	%
	% Only X and Y are mandatory arguments
	% useLogScale = 1 to take the log of the data (to visalise data where
    %   there are large variations in density)
	% axisLimits = [x_min, x_max, y_min, y_max], [] to omit
	% bins = [number_x_bins, number_y_bins] (defines resolution of output)
    %
    % e.g. cloudPlot(x,y,[],1,[500 500]) to produce a cloudPlot with 
    %       500x500 bins, displaying the log of the density with automatic
    %       axis limits
	%
	% By Romesh Abeysuriya, based heavily on cloudPlot.m by Daniel Armyr
    %
    % (from Daniel Armyr's cloudPlot.m):
    % A cloudplot is in essence a 2 dimensional histogram showing the 
    % density distribution of the data described by X and Y. 
    % As the plot displays density information, the
    % dimensionality of X and Y are ignored. The only requirement is that X and
    % Y have the same number of elements. Cloudplot is written to visualize
    % large quantities of data and is most appropriate for datasets of 10000
    % elements or more.
    % suppress = 1 will not render the plot
    % bw = 1 will do canvas = canvas>0 so you get a dot or not (no density)
    
    if strcmp(class(X),'struct') % If we are simply showing an existing plot
        plotdata = X;
        show_canvas(plotdata);
        return
    end

	% Check the data size
	assert (numel(X) == numel(Y),'The number of elements in X and Y must be the same');
	if nargin < 8 || isempty(bw); bw = 0; end
	if nargin < 7 || isempty(suppress); suppress = 0; end
	if nargin < 6 || isempty(colour); colour = []; end
	if nargin < 4 || isempty(useLogScale); useLogScale = false; end
    if nargin < 3 || isempty(axisLimits) % Calculate default axis limits
		 axisLimits = [min(X(:)) max(X(:)), min(Y(:)) max(Y(:))] ;
    end
    
	%Remove any nans or Infs in the data
	pointSelect = isfinite(X) & isfinite(Y) & X >= axisLimits(1) & X <= axisLimits(2) & Y >= axisLimits(3) & Y <= axisLimits(4);
	X = X(pointSelect);
	Y = Y(pointSelect);
	
	set(gca,'Xlim',[axisLimits(1) axisLimits(2)],'Ylim',[axisLimits(3) axisLimits(4)],'units','pixels');
	
    if nargin < 5 || isempty(bins) % Calculate the default number of pixels
		plotsize = get(gca,'Position');
		bins = ceil(plotsize(3:4));
	end
	binSize(2) = diff(axisLimits(3:4))./(bins(2)-1);
	binSize(1) = diff(axisLimits(1:2))./(bins(1)-1);
	canvas = zeros(bins);
    if ~isempty(colour)
        canvas2 = canvas;
    end
    
	% Calculate bin indices for each of the points in the data arrays
	xBinIndex = floor((X - axisLimits(1))/binSize(1))+1; 
	yBinIndex = floor((Y - axisLimits(3))/binSize(2))+1;

	% Collect the indices
	for i = 1:numel(xBinIndex);
		canvas(xBinIndex(i),yBinIndex(i)) = canvas(xBinIndex(i),yBinIndex(i)) + 1;
		if ~isempty(colour) % If a colour axis was provided
		    canvas2(xBinIndex(i),yBinIndex(i)) = canvas(xBinIndex(i),yBinIndex(i)) + colour(i);
		end
	end

	if ~isempty(colour)  % Canvas is mean values
	    canvas = canvas2./canvas;
	else % Do some postprocessing
		% Warning- if canvas>0 is used, log10 must not be used as log10(1)=0
	    if useLogScale
		    canvas = log10(canvas);
        end
     	if bw
        	canvas = canvas > 0; % Show binary rather than graded data
        end
	    colormap gray % Use grayscale so background is black
	    canvas  = imcomplement(canvas); % Display a negative image
	end
	
	plotdata.canvas = canvas;
    plotdata.axisLimits = axisLimits;
    plotdata.colormap = colormap;
    
    if suppress == 0
        show_canvas(plotdata);
    end

end

function show_canvas(plotdata)
	canvas= plotdata.canvas;
    axisLimits = plotdata.axisLimits;
    if isempty(canvas)
        return
    end
    colormap(plotdata.colormap);
	set(gca,'Xlim',[axisLimits(1) axisLimits(2)],'Ylim',[axisLimits(3) axisLimits(4)],'units','pixels');
	image_handle = imagesc([axisLimits(1) axisLimits(2)],[axisLimits(3) axisLimits(4)],canvas'); % Show the image
	axis('xy') % Reverse the y-axis
	set(gca,'units','normalized' ); % Enable resizing
end
	





