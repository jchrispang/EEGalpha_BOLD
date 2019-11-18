function widget(widget_mode,initial_pars,input_p, plot_ss)
	% This program shows the graphical interface widget
	% The parameters are stored internally in the figure guidata
	% To retrieve the params object corresponding to the current contents
	% of the window, use
	% 	h = guidata(gcf);
	% 	p = h.p;
	% 
	% USAGE
	% widget_mode - Optionally chooose from 'nus','gab','xyz'
	% initial_pars - Optionally specify the initial values for the sliders
	% input_p - Provide a model.params object to change qmax, theta etc. 
    % plot_ss -- falg to plot steady state equation instead of erp
    % model.widget('nus',[],ec, 1) - to plot steady states

    
	if nargin < 3 || isempty(input_p)
		input_p = model.params;
        h.plot_ss = 0;
	end

	if nargin < 1 || isempty(widget_mode)
		h.widget_mode = 'xyz';
        h.plot_ss = 0;
	else
		h.widget_mode = widget_mode;
    end
    
    if ~strcmp(h.widget_mode, 'nus')
        h.plot_ss = 0;
    else
        h.plot_ss = plot_ss;
    end

	switch h.widget_mode
		case 'nus'
			h.param_symbols = {'\nu_{ee}','\nu_{ei}','\nu_{es}','\nu_{se}','\nu_{sr}','\nu_{sn}','\nu_{re}','\nu_{rs}','\alpha','\beta','t_0'};
			h.param_units = {'','','','','','ms','ms','ms'};
			h.limits = [0.3    -60    0.02  0.03  -32     0.01  0.07  0.02     10      100    0.075 ;...
		   		        60   -0.04   40    48    -0.02   40    6    7	      100    800    0.14    ];
		   	h.initial_pars = [1.5254   -3.0228    0.5675    3.4474   -1.4651    3.5933    0.1696    0.0507  83.3333  769.2308 0.0850];
		case 'gab'
			h.param_symbols = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta','t_0'};
			h.param_units = {'','','','','','ms','ms','ms'};
			h.limits = [ eps  -40        eps     -40      -14      10      100    0.075 ;...
		   		        20   -eps       40      -eps      -eps      100    800    0.14    ];
		   	h.initial_pars = [  2.0743   -4.1104    5.9943   -1.6712   -0.6474 83.3333  769.2308 0.0850];

		case 'xyz'
			h.param_symbols = {'X','Y','Z','\alpha','\beta','t_0'};
			h.param_units = {'','','','ms','ms','ms'};
			h.limits = [ eps  -1   eps       10      100    0.075 ;...
		   		        1     1     1      100      800    0.14   ];
		   	h.initial_pars = [0.4059    0.5135    0.0571 83.3333  769.2308 0.0850];
	end

	if nargin >= 2 && ~isempty(initial_pars)
		if length(h.initial_pars) ~= length(initial_pars)
			error('Initial parameters do not have the correct size')
		else
			h.initial_pars = initial_pars;
		end
	end

	h = make_components(h,input_p, h.plot_ss);
	guidata(h.fig,h)
	redraw_callback(h.fig)

function redraw_callback(fig)
	h = guidata(fig);
	for j = 1:length(h.param_symbols)
		set(h.slider_label(j),'String',sprintf('%s = %.4f',h.param_symbols{j},get(h.slider(j),'Value')));
		pars(j) = get(h.slider(j),'Value');
	end

	% Put the parameters into the params object
	switch h.widget_mode
		case 'nus'
			h.p.alpha(1:end) = pars(9);
			h.p.beta(1:end) = pars(10);
			h.p.t0 = pars(11);
			h.p.taues = h.p.t0/2;
			h.p.tause = h.p.t0/2;
			h.p.nus = pars(1:8)/1000;
            
            % Keep all steady states and their stability
            h.ve_root.ve_limits = [-10.0, 10.4];
            h.p.complete_gab([],h.ve_root.ve_limits,linspace(0,200,100))
            h.ve_root.phie = nan(5, 1);
            h.ve_root.steady_states = nan(5,1);
            h.ve_root.ves = nan(5,1);
            h.ve_root.stability = -ones(5, 1);
            
            if length(h.p.stab) > 1, % more than one solution
                h.ve_root.phie = h.p.phia(:, 1);
                h.ve_root.steady_states(1:h.p.nroots) = 0;
                h.ve_root.stability(1:h.p.nroots) = h.p.stab;
                h.ve_root.ves(1:h.p.nroots) = h.p.ve;
            elseif isempty(h.p.stab) % no solution -- something went wrong
                disp('Did not find steady states - something went wrong')
            elseif h.p.stab == 1,% one solution
                h.ve_root.phie(1) = h.p.phia(1);
                h.ve_root.steady_states(1) = 0;
                h.ve_root.ves(1) = h.p.ve;
                h.ve_root.stability(1) = h.p.stab;
            end
            
            % If there are multiple roots, this bit of code gets rid of
            % them but the first.
			if length(h.p.stab) > 1 % On the first run, or on returning to stability, we could find more than one solution
				% Pick the first stable solution - this will lead to the widget missing cases where there is more than one solution
                idx = find(h.p.stab,1,'first');
				if isempty(idx)
					idx = 1;
				end
				h.p.phia = h.p.phia(idx,:);
				h.p.complete_gab(1);
            end
            
			[stab,f,P,stab_w,q2] = h.p.dispersion_calc(1);
			if ~isfinite(stab)
				h.p.xyz = [NaN NaN NaN];
				erp_t = NaN;
				erp_V = NaN;
				set(h.q2.steady_state_text,'Visible','on');
			else
				set(h.q2.steady_state_text,'Visible','off');
				[erp_t,erp_V] = h.p.erp;
            end

            
            if h.plot_ss,
              ve  = linspace(h.ve_root.ve_limits(1),h.ve_root.ve_limits(2),1e5);
              gve = model.ve_root(ve,h.p.nus, h.p.theta, h.p.sigma, h.p.qmax);
            end
		case 'gab'
			h.p.alpha(1:end) = pars(6);
			h.p.beta(1:end) = pars(7);
			h.p.t0 = pars(8);
			h.p.taues = h.p.t0/2;
			h.p.tause = h.p.t0/2;
			h.p.gabcd = pars(1:5);
			h.p.gab = []; % Ensure that the gab are not used
			h.p.xyz = get_xyz(h.p);
			[stab,f,P,stab_w,q2] = h.p.dispersion_calc(1);
			[erp_t,erp_V] = h.p.erp;
		case 'xyz'
			h.p.alpha(1:end) = pars(4);
			h.p.beta(1:end) = pars(5);
			h.p.t0 = pars(6);
			h.p.taues = h.p.t0/2;
			h.p.tause = h.p.t0/2;
			h.p.xyz = pars(1:3);
			[stab,f,P,stab_w,q2] = h.p.xyz_spectrum;
			[erp_t,erp_V] = h.p.xyz_erp;
	end
	% First, compute some stats
	%a = model.get_pkf(f,P);
	%beta_f = a.beta_maxf/a.alpha_maxf
	alpha_freq = 1/(h.p.t0 + 1/h.p.alpha(1) + 1/h.p.beta(1)); % Eq 21 HBM 2004
	sigma_freq = sqrt(h.p.beta(1)*h.p.alpha(1))/2/pi; % Robinson 2002

	% DRAW SPECTRUM
	set(h.spec.spectrum,'XData',f,'YData',P);
	yl = get(h.ax_spec,'YLim');
	set(h.spec.t0_marker(1),'XData',[alpha_freq,alpha_freq],'YData',yl);
	set(h.spec.t0_marker(2),'XData',[2*alpha_freq,2*alpha_freq],'YData',yl);
	set(h.spec.t0_marker(3),'XData',[3*alpha_freq,3*alpha_freq],'YData',yl);
	set(h.spec.spindle_marker,'XData',[sigma_freq,sigma_freq],'YData',yl);
	set(h.spec.spindle_text,'Position',[sigma_freq,yl(1) 1]);
	set(h.spec.t0_text(1),'Position',[alpha_freq,yl(1) 1]);
	set(h.spec.t0_text(2),'Position',[2*alpha_freq,yl(1) 1]);
	set(h.spec.t0_text(3),'Position',[3*alpha_freq,yl(1) 1]);

	% DRAW TENT
	try
		[x,y,z,u] = tent.compute_mex(h.p.alpha(1),h.p.beta(1),h.p.t0,h.p.gammae);
	catch
		[x,y,z,u] = tent.compute(h.p.alpha(1),h.p.beta(1),h.p.t0,h.p.gammae);
	end

	set(h.tent.surf,'XData',x,'YData',y,'ZData',z,'CData',u);
	set(h.tent.surf,'LineWidth',0.5,'FaceColor','none')
	edge_alphadata = x+y<0.90 & x > 0 & x < 1 & y > -1;
	set(h.tent.surf,'AlphaData',+edge_alphadata,'EdgeAlpha','flat','AlphaDataMapping','none')

	set(h.tent.marker,'XData',h.p.xyz(1),'YData',h.p.xyz(2),'ZData',h.p.xyz(3));
	set(h.tent.drop,'XData',[h.p.xyz(1) h.p.xyz(1)],'YData',[h.p.xyz(2) h.p.xyz(2)],'ZData',[0 h.p.xyz(3)]);
	set(h.tent.a0a1,'XData',[0 1+2/(h.p.gammae*h.p.t0) NaN 0 1+2/(h.p.gammae*h.p.t0)],'YData',[-2/(h.p.gammae*h.p.t0) -2/(h.p.gammae*h.p.t0) NaN 1 -2/(h.p.gammae*h.p.t0)],'ZData',[0 0 NaN 0 0]);

	% DRAW ERP or DRAW VE_ROOT
    if h.plot_ss > 0,
        set(h.ve_root.line, 'XData', ve, 'YData', gve)
            for ii=1:length(h.ve_root.ves)
               if h.ve_root.stability(ii)==1,
                   set(h.ve_root.ss_marker(ii), 'XData', h.ve_root.ves(ii),  'YData', h.ve_root.steady_states(ii), 'markersize', 7, 'linewidth', 2, 'Color', [0 0 1]);
               elseif h.ve_root.stability(ii)==0,
                   set(h.ve_root.ss_marker(ii), 'XData', h.ve_root.ves(ii),  'YData', h.ve_root.steady_states(ii), 'markersize', 7, 'linewidth', 2, 'Color', [1 0 0]);
               else
                   set(h.ve_root.ss_marker(ii), 'XData', h.ve_root.ves(ii),  'YData', h.ve_root.steady_states(ii), 'markersize', 7, 'linewidth', 2);
               end
            end
        set(h.ve_root.root_marker, 'XData', [h.p.ve, h.p.ve], 'YData', [-500, 500])
        h.erp.line = [];
    else
        set(h.erp.line,'XData',erp_t,'YData',erp_V)
    end
	% DRAW Q2
	set(h.q2.line,'XData',real(q2),'YData',imag(q2));

	% Set stability colours
	if isfinite(stab) && stab
		set([h.spec.spectrum h.q2.line h.erp.line],'Color','b')
		set(h.tent.marker,'MarkerFaceColor','g','CData',[0 1 0])
	else
		set([h.spec.spectrum h.q2.line h.erp.line],'Color','r')
		set(h.tent.marker,'MarkerFaceColor','r','CData',[1 0 0])
	end

function xyz = get_xyz(p)
	xyz(1) = p.gabcd(1)./(1-p.gabcd(2));
	xyz(2) = (p.gabcd(3)+p.gabcd(4))./(1-p.gabcd(5))./(1-p.gabcd(2));
	xyz(3) = -p.gabcd(5).*p.alpha(1).*p.beta(1)./((p.alpha(1)+p.beta(1)).^2);

function h = prepare_spectrum(h)
	h.spec.spectrum = loglog(h.ax_spec,NaN,NaN,'b');
	set(h.ax_spec,'XLim',[0.25 45]);
	hold(h.ax_spec,'on');
	xlabel(h.ax_spec,'Frequency (Hz)');
	ylabel(h.ax_spec,'Power (normalized)');

	h.spec.spindle_marker = plot(h.ax_spec,[NaN NaN],[NaN NaN],'m--');
	h.spec.t0_marker(1) = plot(h.ax_spec,[NaN NaN],[NaN NaN],'g--');
	h.spec.t0_marker(2) = plot(h.ax_spec,[NaN NaN],[NaN NaN],'g--');
	h.spec.t0_marker(3) = plot(h.ax_spec,[NaN NaN],[NaN NaN],'g--');

	h.spec.spindle_text = text(NaN,NaN,'\sigma','VerticalAlignment','bottom','Parent',h.ax_spec);
	h.spec.t0_text(1) = text(NaN,NaN,'f_\alpha','VerticalAlignment','bottom','Parent',h.ax_spec);
	h.spec.t0_text(2) = text(NaN,NaN,'2f_\alpha','VerticalAlignment','bottom','Parent',h.ax_spec);
	h.spec.t0_text(3) = text(NaN,NaN,'3f_\alpha','VerticalAlignment','bottom','Parent',h.ax_spec);

	hold(h.ax_spec,'off');

function h = prepare_tent(h)
	try
		[x,y,z,u] = tent.compute_mex(h.p.alpha(1),h.p.beta(1),h.p.t0,h.p.gammae);
	catch
		[x,y,z,u] = tent.compute;
		fprintf(2,'Compile tent.compute_mex for faster performance\n');
	end

	h.tent.surf = mesh(h.ax_tent,x,y,z,u);
	hold(h.ax_tent,'on')
	h.tent.marker = scatter3(h.ax_tent,NaN,NaN,NaN,50,'go','MarkerFaceColor','g');
	h.tent.drop = plot3(h.ax_tent,[NaN NaN],[NaN NaN],[NaN NaN],'r--');
	h.tent.a0a1 = plot3(h.ax_tent,NaN,NaN,NaN,'k--','LineWidth',2);



	colorbar('peer',h.ax_tent,'Location','southoutside')
	set(h.ax_tent,'CLim',[0 45],'ZLim',[0 1],'YLim',[-1 1],'View',[100  20],'XLim',[min(x(:)) max(x(:))],'YLim',[min(y(:)) max(y(:))],'ZLim',[0 1])
	xlabel(h.ax_tent,'X');
	ylabel(h.ax_tent,'Y');
	zlabel(h.ax_tent,'Z');
	hold(h.ax_tent,'off')

function h = prepare_q2re2(h)
	h.q2.line = plot(h.ax_q2,NaN,NaN);
	hold(h.ax_q2,'on')
	plot(h.ax_q2,[-800 400],[0,0],'k--','XLimInclude','off')
	plot(h.ax_q2,[0,0],[-800 400],'k--','YLimInclude','off')
	xlabel(h.ax_q2,'Re(q^2)');
	ylabel(h.ax_q2,'Im(q^2)');
	set(h.ax_q2,'XLim',[-800 400],'YLim',[-800 400])
	h.q2.steady_state_text = text(0,0,'NO STEADY STATES FOUND','FontWeight','bold','Color','r','HorizontalAlignment','right','Parent',h.ax_q2,'VerticalAlignment','bottom','Visible','off');
	hold(h.ax_q2,'off')


function h = prepare_erp(h)
	h.erp.line = plot(h.ax_erp,NaN,NaN);
	set(h.ax_erp,'XLim',[0 1],'YLim',[-10 50]);
	xlabel(h.ax_erp,'Time');
	ylabel(h.ax_erp,'Voltage');
    
function h = prepare_ve_root(h)
    h.ve_root.line = plot(h.ax_ve_root,NaN,NaN);
    h.ve_root.steady_states = NaN;
    h.ve_root.ves  = NaN;
    h.ve_root.phie = NaN; 
    h.ve_root.stability = nan(5, 1);
    hold(h.ax_ve_root, 'on')
    h.ve_root.ss_marker(1) = plot(h.ax_ve_root, h.ve_root.ves,  h.ve_root.steady_states, 'ko');
    h.ve_root.ss_marker(2) = plot(h.ax_ve_root, h.ve_root.ves,  h.ve_root.steady_states, 'ko');
    h.ve_root.ss_marker(3) = plot(h.ax_ve_root, h.ve_root.ves,  h.ve_root.steady_states, 'ko');
    h.ve_root.ss_marker(4) = plot(h.ax_ve_root, h.ve_root.ves,  h.ve_root.steady_states, 'ko');
    h.ve_root.ss_marker(5) = plot(h.ax_ve_root, h.ve_root.ves,  h.ve_root.steady_states, 'ko');
    h.ve_root.root_marker  = plot(h.ax_ve_root,[NaN NaN],[NaN NaN],'g--');
    plot(h.ax_ve_root, [-2 2], [0 0], 'k--','XLimInclude','off')
    set(h.ax_ve_root, 'YLim', [-100, 200], 'XLim', [-0.04, 0.04])
    xlabel(h.ax_ve_root, 've')
    ylabel(h.ax_ve_root, 'g(ve)')
    hold(h.ax_ve_root,'off')

        
function h = make_components(h,input_p, plot_ss)
    % plot_ss = 1, plots the steady state equation instead of erp
    if nargin < 3,
        plot_ss = 0;
    end
	h.fig = figure;
	h.hg2 = isa(h.fig,'handle');
	if ~h.hg2
		fprintf(2,'Use Matlab R2014b or newer for optimal appearance\n');
	end

	pos = get(h.fig,'Position');
	set(h.fig,'Position',pos.*[1 1 2 1.5]-pos.*[0.5 0.5 0 0]); % Twice as wide, and a bit bigger too
	set(h.fig,'Color',[0.94 0.94 0.94])

	h.ax_tent = axes('Parent',h.fig,'Position',[0.06 0.58 0.3 0.4]); 
	h.ax_q2 = axes('Parent',h.fig,'Position',[0.06 0.08  0.3 0.4]); 
	h.ax_spec = axes('Parent',h.fig,'Position',[0.45 0.58  0.3 0.4]);
    if plot_ss,
        h.ax_ve_root = axes('Parent',h.fig,'Position',[0.45 0.08  0.3 0.4]); 
    else
	    h.ax_erp = axes('Parent',h.fig,'Position',[0.45 0.08  0.3 0.4]); 
    end
	h.panel_sliders = uipanel('Parent',h.fig,'Position',[0.79 0.05 0.2 0.9],'Title','PARAMETERS'); 
	
	h.p = input_p;
    h.plot_ss = plot_ss;

	slider_height = min(0.85/length(h.param_symbols)/2,0.05);
	gap = 0.05;
	for j = 1:length(h.param_symbols)
		if h.hg2
			h.slider_label(j) = annotation(h.panel_sliders,'textbox','EdgeColor','none','HorizontalAlignment','center','Units','normalized','String','INIT','Position',[0.05 1-gap*j-slider_height*(j-1)-0.05  0.9 slider_height],'Interpreter','tex');
		else
			h.slider_label(j) = uicontrol('Style','text','Parent',h.panel_sliders,'HorizontalAlignment','center','Units','normalized','String','INIT','Position',[0.05 1-gap*j-slider_height*(j-1)-0.05  0.9 slider_height]);
		end
		h.slider(j) = uicontrol('Style','slider','Units','normalized','Position',[0.05 1-gap*j-slider_height*(j-1)  0.9 slider_height],'Parent',h.panel_sliders,'Max',h.limits(2,j),'Min',h.limits(1,j),'Value',h.initial_pars(j));
		h.listener(j) = addlistener(h.slider(j),'Value','PostSet',@(s,e)  redraw_callback(h.fig) );
	end

	% Now prepare each of the components
	h = prepare_spectrum(h);
	h = prepare_tent(h);
	h = prepare_q2re2(h);
    
    if plot_ss,
        h = prepare_ve_root(h);
    else
        h = prepare_erp(h);
    end
	colormap([0.2081    0.1663    0.5292 % Copied from Matlab 2015a)
    0.2116    0.1898    0.5777
    0.2123    0.2138    0.6270
    0.2081    0.2386    0.6771
    0.1959    0.2645    0.7279
    0.1707    0.2919    0.7792
    0.1253    0.3242    0.8303
    0.0591    0.3598    0.8683
    0.0117    0.3875    0.8820
    0.0060    0.4086    0.8828
    0.0165    0.4266    0.8786
    0.0329    0.4430    0.8720
    0.0498    0.4586    0.8641
    0.0629    0.4737    0.8554
    0.0723    0.4887    0.8467
    0.0779    0.5040    0.8384
    0.0793    0.5200    0.8312
    0.0749    0.5375    0.8263
    0.0641    0.5570    0.8240
    0.0488    0.5772    0.8228
    0.0343    0.5966    0.8199
    0.0265    0.6137    0.8135
    0.0239    0.6287    0.8038
    0.0231    0.6418    0.7913
    0.0228    0.6535    0.7768
    0.0267    0.6642    0.7607
    0.0384    0.6743    0.7436
    0.0590    0.6838    0.7254
    0.0843    0.6928    0.7062
    0.1133    0.7015    0.6859
    0.1453    0.7098    0.6646
    0.1801    0.7177    0.6424
    0.2178    0.7250    0.6193
    0.2586    0.7317    0.5954
    0.3022    0.7376    0.5712
    0.3482    0.7424    0.5473
    0.3953    0.7459    0.5244
    0.4420    0.7481    0.5033
    0.4871    0.7491    0.4840
    0.5300    0.7491    0.4661
    0.5709    0.7485    0.4494
    0.6099    0.7473    0.4337
    0.6473    0.7456    0.4188
    0.6834    0.7435    0.4044
    0.7184    0.7411    0.3905
    0.7525    0.7384    0.3768
    0.7858    0.7356    0.3633
    0.8185    0.7327    0.3498
    0.8507    0.7299    0.3360
    0.8824    0.7274    0.3217
    0.9139    0.7258    0.3063
    0.9450    0.7261    0.2886
    0.9739    0.7314    0.2666
    0.9938    0.7455    0.2403
    0.9990    0.7653    0.2164
    0.9955    0.7861    0.1967
    0.9880    0.8066    0.1794
    0.9789    0.8271    0.1633
    0.9697    0.8481    0.1475
    0.9626    0.8705    0.1309
    0.9589    0.8949    0.1132
    0.9598    0.9218    0.0948
    0.9661    0.9514    0.0755
    0.9763    0.9831    0.0538])
