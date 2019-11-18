function [bins,binMax] = bin_max(x,bins)
	% Find the maximum value of binned data
	%
	% bin_max(x,nbins)
	%
	% x is a two column matrix
	% For bins of the first column, get the max value of the second
	% This is good for finding boundaries
	% Much code taken from 
	% http://blogs.mathworks.com/videos/2009/01/07/binning-data-in-matlab/
	if numel(bins) > 1
		binEdges = bins;
	else
		binEdges = linspace(min(x(:,1)),max(x(:,1)),bins+1);
        bins = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
	end

	nbins = length(binEdges);
	binMax = zeros(1,nbins);

	[h,whichBin] = histc(x(:,1),binEdges);
	for i = 1:nbins
		flagBinMembers = (whichBin == i);
		if ~any(flagBinMembers)
            binMax(i) = 0;
        else
            binMax(i) = max(x(flagBinMembers,2));
        end
	end
	
	
	%scatter(bins,binMax,'bx');
	
	
