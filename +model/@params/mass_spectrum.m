function varargout = mass_spectrum(p,varargin)
	if (nargin < 2 || isempty(verbose)) && nargout == 0
	    verbose = 1; % If user requested no outputs, draw a plot
	elseif nargin < 2 || isempty(verbose)
	    verbose = 0;
	end

    [f,P] = p.spectrum(verbose,[],[],true);

    if nargout > 0
    	varargout{1} = f;
    	varargout{2} = P;
    end
