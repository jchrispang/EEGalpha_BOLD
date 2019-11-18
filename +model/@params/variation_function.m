function fcn = variation_function(self,var_type,varargin)
	% This function serves as a container for various patterns of spatial make_variation
	% that can be generically applied to different parameters by wrapper functions
	% varargin{1} is always the current value of the parameter in question
	% then the remaining arguments were provided by the parent or user
	switch var_type
		case 'cosine'
			fcn = cosine(self.Ly,varargin{1:end}); 
		case 'cosine_phased'
			fcn = cosine_phased(self.Ly,varargin{1:end}); 
		case 'cosine_2d_phased'
			fcn = cosine_2d_phased(self.Lx,self.Ly,varargin{1:end}); 
		case 'gaussian'
			fcn = gaussian(varargin{1:end});
		case 'single_step'
			fcn = single_step(self.Ly,varargin{1:end});
		otherwise
			error('Variation type not found')
	end
end

function fcn = cosine(Ly,par,val)
	% One variation across the length of the scalp
	% The inverted sign means that the parameter is MAXIMUM at the middle
	fcn = @(x,y) par + val*cos((y-Ly/2)*2*pi/Ly);
end

function fcn = cosine_phased(Ly,par,val,ymax)
	% One variation across the length of the scalp
	% Inverted sign means that the parameter is maximum at position ymax 
	fcn = @(x,y) par + val*cos((y-ymax)*2*pi/Ly);
end

function fcn = cosine_2d_phased(Lx,Ly,par,x_amp,y_amp,xmax,ymax)
	% One variation across the length of the scalp
	% Inverted sign means that the parameter is maximum at position ymax 
	% Don't forget that the amp parameters need to be allowed to become negative if required
	fcn = @(x,y) par + x_amp*cos((x-xmax)*2*pi/Lx) + y_amp*cos((y-ymax)*2*pi/Ly);
end

function fcn = gaussian(par,p1,x1,y1,sig)
	fcn = @(x,y) par + p1.*exp(( -(x-x1).^2 -(y-y1).^2)/(2*sig.^2));
end

function fcn = single_step(Ly,par,offset,yval)
	fcn = @(x,y) par + offset*(y>yval);
end