function varargout = bin_cloud(state,colour,colour_param_string)
    % Make the ubiquitous cloudplot
    % colour_param_string is something like data.gab_final(:,8) or something that
    % evaluates to a vector of values the same size as the data. data.spec_analysis.int_alpha for example

    if nargin < 2 || isempty(colour)
        colour = [];
    end

    if nargin < 1 || isempty(state)
        error('Need input state or data');
    end

    if isstruct(state)
        data = state;
    else
        switch state
            case {'ec','eo'}
                a = load('pdb_wake');
                b = load('pdb_wake_analysis_new');
            case {'rem','n1','n2','n2s','n3'}
                a = load('pdb_sleep');
                b = load('pdb_sleep_analysis_new');
        end

        data = catstruct(a,b);
        c = load('pdb_specvalid');
        data.alllim = c.alllim.(state);
    end

    if nargin < 3 || isempty(colour_param_string)
        colour_param = [];
    else
        colour_param = eval(colour_param_string);
    end

    [alllim,bin_3d_data] = bin_3d(data,colour_param);
    
    getbinspectra_data.gab = data.gab_final(alllim,:);
    getbinspectra_data.xyz = data.xyz_final(alllim,:);
    getbinspectra_data.phi = data.phi_final(alllim,:);
    getbinspectra_data.nus = data.nus_final(alllim,:);
    getbinspectra_data.ve = data.ve_final(alllim,:);
    getbinspectra_data.state = state;

    if ~isempty(colour) % A colour was specified
        set(gca,'CLim',[0 1])
        switch colour
            case 'r'
                cval = 0.9;
            case 'b'
                cval = 0.1;
            case 'c'
                cval = 0.38;
            case 'g'
                cval = 0.5;
            case 'y'
                cval = 0.62;
            case 'o'
                cval = 0.7;
            case 'm'
                cval = 0.3;
        end
        bin_3d_data.cdata(bin_3d_data.alphadata>0) = cval; % Make all points the same colour
    end

    vol = vol3d_alpha(bin_3d_data.xyzlim,bin_3d_data.cdata,+(bin_3d_data.alphadata>0));

    set(gca,'ALim',[0 4]);

    u.plot_type = 'bin_cloud';
    h = get(gca,'Children');
    h = get(h(end),'CData');
    u.xl = get(gca,'XLim');
    u.yl = get(gca,'YLim');
    u.zl = get(gca,'ZLim');
    u.gridres = diff(u.xl)/size(h,2);
    u.xbins = u.xl(1):u.gridres:u.xl(2);
    u.ybins = u.yl(1):u.gridres:u.yl(2);
    u.zbins = u.zl(1):u.gridres:u.zl(2);
    u.view = get(gca,'View');
    u.surface_handles = vol.surface_handles;

    set(gcf,'Userdata',u);

    if nargin < 2 || isempty(colour)% Only draw a colourbar if it is meaningful
        colorbarhandle = colorbar;
    end
    
    if nargout > 0
        varargout{1} = getbinspectra_data;
        varargout{2} = vol;
    end
