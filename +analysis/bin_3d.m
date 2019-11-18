function varargout = bin_3d(data,colour_parameter)
    % outputs: [alllim,bin_3d_data] = bin_3d(data,colour_parameter)

    % Bin the data
    xyzlim = [0  -1     0 ;...
	          1   1   1.3];

    gridres = 0.025;

    xlim = (data.xyz_final(:,1) < xyzlim(2,1) & data.xyz_final(:,1) > xyzlim(1,1));
    ylim = (data.xyz_final(:,2) < xyzlim(2,2) & data.xyz_final(:,2) > xyzlim(1,2));
    zlim = (data.xyz_final(:,3) < xyzlim(2,3) & data.xyz_final(:,3) > xyzlim(1,3));
    gratio = data.gab_final(:,1)./(data.gab_final(:,2));
    
    if ~isfield(data,'alllim')
        data.alllim = (xlim & ylim & zlim) & data.stab_final==1 & gratio < -1/2 & data.specvalid==0;
    end

    fprintf('Displaying %d points\n',sum(data.alllim));

    xyz = data.xyz_final(data.alllim,:); 

    if nargin < 2 || isempty(colour_parameter)
        % Choose which parameter will be used for colour
        %olour_parameter = ones(size(data.xyz_final(:,1)));
        %colour_parameter = stab_final(data.alllim); % stability
        %colour_parameter = phi_final(data.alllim,1); % one of the firing rates
        %colour_parameter = gratio(data.alllim); % the gee/gei ratio
        colour_parameter = data.gab_final(data.alllim,1); % one of the nuab
        %colour_parameter = sum(nus_final(data.alllim,1:2),2);
        %colour_parameter = gab_final(data.alllim,1); % one of the gab
        %colour_parameter = pkf_final(data.alllim,15); % one of the pkf parameters
        %colour_parameter = maxfreq(data.alllim); % frequency of the strongest peak
        %colour_parameter = phi_final(data.alllim,3).*nus_final(data.alllim,3);
        %colour_parameter = sum(nus_final(data.alllim,1:2),2);
        %colour_parameter = gab_final(data.alllim,5).*gab_final(data.alllim,7) + gab_final(data.alllim,4);
        %colour_parameter = spec_analysis.hf_slope(data.alllim);
        %colour_parameter = 1 - xyz_final(data.alllim,1) - xyz_final(data.alllim,2);
        %colour_parameter = data.spec_analysis.lf_slope(data.alllim); % Total sigma power
        %colour_parameter = gab_final(data.alllim,1)+gab_final(data.alllim,2);
        %colour_parameter = validate_problems(:,10);
    else
        colour_parameter = colour_parameter(data.alllim);
    end
    % There used to be some other code to get (xyz)BinIndex, but it hasn't 
    % been used for a while, so I have removed it
    xBinIndex = floor((xyz(:,1)-xyzlim(1,1))./gridres)+1;
    yBinIndex = floor((xyz(:,2)-xyzlim(1,2))./gridres)+1;
    zBinIndex = floor((xyz(:,3)-xyzlim(1,3))./gridres)+1;

    V = zeros((xyzlim(2,2)-xyzlim(1,2))/gridres,(xyzlim(2,1)-xyzlim(1,1))/gridres,(xyzlim(2,3)-xyzlim(1,3))/gridres);
    Vcount = V;

    for zval = 1:max(zBinIndex)
	    this_z = zBinIndex==zval;
	
	    xb = xBinIndex(this_z);
	    yb = yBinIndex(this_z);
        colour_z = colour_parameter(this_z);

	    for j = 1:length(xb)-1
	        %V(yb(j),xb(j),zval) =  max(V(yb(j),xb(j),zval),colour_z(j)); 

            % Colour using the to_add parameter. If it is not finite, preserve the average
            if ~isfinite(colour_z(j))
                V(yb(j),xb(j),zval) =  V(yb(j),xb(j),zval) * Vcount(yb(j),xb(j),zval) / ( Vcount(yb(j),xb(j),zval)+1 );
		    else
		        V(yb(j),xb(j),zval) =  V(yb(j),xb(j),zval) + colour_z(j); 
	        end
	        
		    % Count (always enable for density transparency)
		    Vcount(yb(j),xb(j),zval) =  Vcount(yb(j),xb(j),zval)+1;
	    end
    end

    V = V./Vcount; % Enable with finding mean properties

    bin_3d_data.xyzlim = xyzlim;
    bin_3d_data.gridres = gridres;
    bin_3d_data.cdata = V;
    bin_3d_data.alphadata = Vcount;

    if nargout >= 1
        varargout{1} = data.alllim;
    end
    if nargout >= 2
        varargout{2} = bin_3d_data;
    end
