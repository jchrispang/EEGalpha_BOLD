function getbinspectra(d)
    % Task argument is used to determine which alpha and beta values
	% Automatically show spectra in a bin
    
    % d is supposed to be the output from bin_cloud
    % at the moment, ugly hack allows d to be a string state
    
    defaults = get(gcf,'UserData');
    if isempty(defaults) || ~(strcmp(defaults.plot_type,'bin_cloud')) % Make a plot if not already open
        warning('You need to have a bin_cloud plot open before running getbinspectra()')
        d = bin_cloud(d);
    end
    
	task_str = d.state;

	tasks = {'all','ec','eo','rem','n1','n2','n2s','n3'};
	wake = [1 1 1 0 0 0 0 0];

	load data/pdb_specvalid validate_fcn;
	task = find(strcmp(task_str,tasks));
	wake = wake(task);
	validate_fcn = validate_fcn.(tasks{task});

	fprintf('Checking %s (wake=%d)\n',tasks{task},wake);


	main_figure = gcf;
	main_axes = gca;

	% STAGE 1- GET A LIST OF INDEXES FOR POINTS TO CYCLE SPECTRA THROUGH
	figure(main_figure);
	cleanupObj = onCleanup(@() getbinspectra_cleanup(defaults,main_figure,d));

	% Display slices of the plot
	set(main_axes,'View',[90,90])
	alphamap([0 ones(1,63)])
	%vol3d_alpha(p);
	j = 2;

	% UNCOMMENT TO GET BACK TO ORIGINAL FUNCTIONALITY
	while j <= length(defaults.zbins)
		set(main_axes,'ZLim',[0 defaults.zbins(j)])
		select_z(j,main_figure)
		title(main_axes,sprintf('XY slices- Z = %2.2f',defaults.zbins(j)))
		drawnow
		brk = input('Accept Z value, or enter a Z value: ', 's')
		if isempty(brk)
			j = j+1;
		elseif str2num(brk)
		    [a,j] = min(abs(defaults.zbins-str2num(brk)));
		else
            break
		end
	end

	disp('Select bin')
	coord = ginput(1)

	if isempty(coord) || coord(1) <= defaults.xbins(1) || coord(1) >= defaults.xbins(end) || coord(2) <= defaults.ybins(1) || coord(2) >= defaults.ybins(end)
		disp('You clicked outside the plot!')
		return
    end
    
    xcorner = find(coord(1)>=defaults.xbins,1,'last');
    ycorner = find(coord(2)>=defaults.ybins,1,'last');
    % Now get the bin edges
	xedges = [defaults.xbins(xcorner),defaults.xbins(xcorner+1)];
	yedges = [defaults.ybins(ycorner),defaults.ybins(ycorner+1)];
	zedges = [defaults.zbins(j-1),defaults.zbins(j)];
	fprintf('SELECTED BIN\nLimit\tMin\tMax\n');
	fprintf('X\t%2.3f\t%2.3f\n',xedges(1),xedges(2));
	fprintf('Y\t%2.3f\t%2.3f\n',yedges(1),yedges(2));
	fprintf('Z\t%2.3f\t%2.3f\n',zedges(1),zedges(2));

	xva = (d.xyz(:,1) >= xedges(1) & d.xyz(:,1) <= xedges(2));
	yva = (d.xyz(:,2) >= yedges(1) & d.xyz(:,2) <= yedges(2));
	zva = (d.xyz(:,3) >= zedges(1) & d.xyz(:,3) <= zedges(2));

	v = find(xva & yva & zva);
	% v is a list of all the indexes of points to iterate over
	
	% STAGE 2: CYCLE THROUGH SPECTRA. NEED TO UPDATE CLEANUP FUNCTION
	new_figure = figure;
	new_axes = gca;
	function index = jfunc()
	    index = v(j);
	end
	cleanupObj = onCleanup(@() getbinspectra_cleanup(defaults,main_figure,d,new_figure,wake,@jfunc));

	assignin('base','found_points',v);
	p1 = params(wake);

	types = {'maxf', 'maxp', 'minf', 'minp', 'pkstr'}; 
	bands = {'delta','theta','alpha','sigma','beta'};

	for j = 1:length(v)
	       	p1.gab = d.gab(v(j),:);
			fprintf(2,'\nENTRY %i of %i(index %i)\n',j,length(v),v(j));
			p1.phin = 5e-5;
			[fa,Pa] = analytic_spectrum(p1,0);
			pkf = get_pkf(fa,Pa);

            figure(new_figure)
            %set(new_figure,'Renderer','painters');
			%loglog(fa,Pa);
			plot(log(fa),log(Pa));
			get_spec_analysis(fa,Pa)
			%set(gca,'XLim',[1 45],'XTick',[1 2 4 6 10 15 25 40]);
			%set(gca,'YLim',[1e-12 1e-6],'YTick',[1e-12 1e-10 1e-8 1e-6])
	   		% validate_fcn(fa,Pa)
		    title(sprintf('X: %1.2f Y: %1.2f Z: %1.2f, %3.2f%% complete',d.xyz(v(j),1),d.xyz(v(j),2),d.xyz(v(j),3),100*j/length(v)));
		    drawnow
		    if filter_points(p1,fa,Pa)
		    	pause
		    end

		    clf(new_figure)
	end
end

function accept = filter_points(p,fa,Pa)
	% Use this function to halt the pausing if a condition is met
	accept = 1;
	return

	if trapz(fa,Pa) > 1e-5 && trapz(fa,Pa) < 3e-5 && pkf(11) > 9
		accept = 1;
	else
		accept = 0;
	end
end

function getbinspectra_cleanup(defaults,main_figure,getbinspectra_data,new_figure,wake,jfunc)
	figure(main_figure)
	select_z([],main_figure)
	set(gca,'ZLim',defaults.zl)
	set(gca,'View',defaults.view)
    alphamap('default')
    title('')

    if nargin < 4 % User cancelled before clicking on a point
    	return
    end

    % Otherwise, create a variable for the currently displayed spectrum
	close(new_figure)
    idx = feval(jfunc);

    fprintf(1,'Index was %i!',idx)
    current_point = params(wake);
    current_point.gab = getbinspectra_data.gab(idx,:);
    current_point.phia = getbinspectra_data.phi(idx,:);
    current_point.ve = getbinspectra_data.ve(idx,:);
    current_point.nus = getbinspectra_data.nus(idx,:);
    current_point.xyz = getbinspectra_data.xyz(idx,:);
    assignin('base','current_point',current_point);
end

function select_z(idx,fig)
	h = get(fig,'UserData');
	h = h.surface_handles{3};

	for j = 1:length(h)
		if isempty(idx) % show everything
			set(h,'Visible','on');
		else
			set(h,'Visible','off');
		end
	end

	if ~isempty(idx)
		set(h(idx),'Visible','on');
	end
end