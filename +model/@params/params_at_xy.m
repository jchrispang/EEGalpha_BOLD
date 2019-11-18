function p = params_at_xy(self,x,y)
	% This function returns a non-spatially varying set of parameters
	% with the values of the spatially varying parameters replacing
	% the spatial variation functions, at the specified coordinates
	p = self.copy();
	p.disable_set = true;

	xyz_mode = 3; % By default, copy existing XYZ

    for j = 1:8
        if ~isempty(self.spatial_alpha) && utils.isfunction(self.spatial_alpha{j})
            p.alpha(j) = self.spatial_alpha{j}(x,y);
        end

        if ~isempty(self.spatial_beta) && utils.isfunction(self.spatial_beta{j})
            p.beta(j) = self.spatial_beta{j}(x,y);
        end

        if ~isempty(self.spatial_nus) && utils.isfunction(self.spatial_nus{j})
            p.nus(j) = self.spatial_nus{j}(x,y);
        end

        if ~isempty(self.spatial_gab) && utils.isfunction(self.spatial_gab{j})
            p.gab(j) = self.spatial_gab{j}(x,y);
            xyz_mode = 2; % Reload XYZ from gab
        end
    end

    for j = 1:5
        if ~isempty(self.spatial_gabcd) && utils.isfunction(self.spatial_gabcd{j})
            p.gabcd(j) = self.spatial_gabcd{j}(x,y);
            xyz_mode = 1; % Reload XYZ from gabcd
        end
    end

    for j = 1:3
        if ~isempty(self.spatial_xyz) && utils.isfunction(self.spatial_xyz{j})
            p.xyz(j) = self.spatial_xyz{j}(x,y);
            xyz_mode = 4; % XYZ has just been computed
        end
    end

	if utils.isfunction(self.spatial_t0)
	    p.t0 = self.spatial_t0(x,y);
	    p.taues = p.t0/2;
	    p.tause = p.t0/2;
	end

	if ~isempty(self.spatial_nus)
		p.complete_gab(1); % Complete gab if the nus have changed
		xyz_mode = 4; % XYZ has just been computed
	end


	% Guard against Gab and Gabcd having conflicting XYZ values
	assert(isempty(self.spatial_gabcd) || isempty(self.spatial_gab)) % XYZ could desync if both are defined. They should never both be defined anyway
	
	if xyz_mode == 3 & isempty(self.xyz)
		if ~isempty(self.gab)
			xyz_mode = 2;
		elseif ~isempty(self.gabcd)
			xyz_mode = 1;
		else
			error('No XYZ source could be identified');
		end
	end

	switch xyz_mode
		case 1
			p.xyz(1) = p.gabcd(1)./(1-p.gabcd(2));
			p.xyz(2) = (p.gabcd(3)+p.gabcd(4))./(1-p.gabcd(5))./(1-p.gabcd(2));
			p.xyz(3) = -p.gabcd(5).*p.alpha(1).*p.beta(1)./((p.alpha(1)+p.beta(1)).^2);
		case 2
			p.xyz(1) = p.gab(1)/(1-p.gab(2));
		    p.xyz(2) = (p.gab(3)*p.gab(4)+p.gab(3)*p.gab(5)*p.gab(7))/((1-p.gab(5)*p.gab(8))*(1-p.gab(2)));
		    p.xyz(3) = -p.gab(5)*p.gab(8)*p.alpha(5)*p.beta(5)/(p.alpha(5)+p.beta(5))^2;
		case 3
			p.xyz = self.xyz;
		case 4
			% Do nothing 
	end

	p.disable_set = false;
