function [reconstruct,alllim] = bin_blob(state,colour,alpha,offset)
        
    if nargin < 4
        % There appear to be artifacts related to both EO and EC terminating at the 
        % X+Y=1 boundary. Therefore, one of them should be slightly offset when 
        % using analysis_to_tent
        offset = 0;
    end

    if nargin < 1 || isempty(state)
        error('Need input state or data')
    end
    
    if nargin < 3 || isempty(alpha)
        alpha = 0.5;
    end
    
    if isstruct(state)
        data = state;
        if nargin < 2 || isempty(colour)
            error('Need to specify a colour if using struct argument')
        end
    else
        cdata = [255 0 0; ...% EO bright red
        0 255 145; ...% EC light green
        255 180 0; ...% REM Orange
        0 159 100; ...% N1 dark green
        0 154 255; ...% N2 light blue
        12 23 178; ...% N3 dark blue
        255 255 0; ...% N2S bright yellow
        ];
        states = {'eo','ec','rem','n1','n2','n3','n2s'};
        iswake = [1 1 0 0 0 0 0];

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
        if nargin < 2 || isempty(colour)
            colour = cdata(this_state,:);
        end
    end

    [alllim,bin_3d_data] = analysis.bin_3d(data);

    % Isosurface method in Rev 785
    x = bin_3d_data.xyzlim(1):bin_3d_data.gridres:bin_3d_data.xyzlim(2,1);
    y = bin_3d_data.xyzlim(1,2):bin_3d_data.gridres:bin_3d_data.xyzlim(2,2);
    z = bin_3d_data.xyzlim(1,3):bin_3d_data.gridres:bin_3d_data.xyzlim(2,3);
    [X,Y,Z] = meshgrid(x(1:end-1),y(1:end-1),z(1:end-1));
    plot_V = smooth3(bin_3d_data.alphadata>0);
    a = data.xyz_final(alllim,:);
    a = a+offset; 
    
    if size(a,1) < 4
        reconstruct = [];
        return
    else
        [V,S] = romesh_utils.alphavol(a,Inf);
    end
    
    k = S.bnd;
    %k = convhull(a(:,1),a(:,2),a(:,3));

    map = unique(k(:));
    map(:,2) = 1:length(map); % First col is the index in k, second row is the new index
    reconstruct.vertices = a(map(:,1),:); % A final list of vertices
    reconstruct.faces = arrayfun(@(b) map(find(map(:,1)==b),2),k);
    reconstruct.colour = colour;
    reconstruct.alpha = alpha;

    %[reconstruct.vertices,reconstruct.faces] = triangulation_subdivide(reconstruct.vertices,reconstruct.faces);
    %[reconstruct.vertices,reconstruct.faces] = triangulation_subdivide(reconstruct.vertices,reconstruct.faces);

    p = trisurf(reconstruct.faces,reconstruct.vertices(:,1),reconstruct.vertices(:,2),reconstruct.vertices(:,3));
    set(p,'FaceColor',reconstruct.colour,'EdgeColor','none','FaceAlpha',reconstruct.alpha,'UserData','blob');

    set(gca,'ALim',[0 1]);
    axis(gca,'equal');
    set(gca,'XLim',[bin_3d_data.xyzlim(1,1),bin_3d_data.xyzlim(2,1)],'YLim',[bin_3d_data.xyzlim(1,2),bin_3d_data.xyzlim(2,2)],'ZLim',[bin_3d_data.xyzlim(1,3),bin_3d_data.xyzlim(2,3)]);
    axis(gca,'vis3d');
    lighting gouraud

    %light('Position',[7.70579 7.10542 8.7023]);
    grid on
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    drawnow
end