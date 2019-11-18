function varargout = neurofield(p,file_id,firemode,int_time,grid_edge,fs,waves,ranseed,fprefix)
    % Run NeuroField via a params object
    % - Accepts a point struct as input
    % - file_id optionally specifies the number for neurofield_*.conf/output
    % - firemode is a 3 element vector that equals
    %   - 0 for fully nonlinear
    %   - 1 for linear
    %   - 2 for quadratic
    %   e.g. [0 0 1] will only use a sigmoid population response in relay nuclei
    %   System is fully nonlinear by default 
    % - grid_edge: specify edge size of the grid (in number of nodes). Center node is automatically outputted
    % - grid_output: outputs grid data for just the excitatory neurons, and returns
    %                spectrum with spectral filtering
    %   This program writes the conf file, executes it, and parses the output
    %
    % WARNING- THIS FUNCTION ASSUMES SAME ALPHA AND BETA FOR ALL CONNECTIONS
    % WARNING- THIS FUNCTION ASSUMES 'neurofield' IS ON THE SHELL PATH
    % WARNING- Depends on nf.read, nf.extract, and nf.grid 
    %
    % Defaults are- 20x20 grid for 0.5m x 0.5m cortex
    %               for 30s at 10000Hz with nonlinearity
    % 
    % PLOTTING BEHAVIOUR
    % If no output arguments are specified >>nf.eirs() then nothing will be
    % produced. The end result is to produce a NeuroField output file
    % If nf = nf.eirs() then nf.read will be called but no plot produced
    % If [nf,f,P] = nf.eirs() then nf.spatial_spectrum will also be called
    % and a plot produced

    if nargin < 9 || isempty(fprefix)
        fprefix = 'neurofield';
    end

    if nargin < 8 || isempty(ranseed)
        ranseed = [];
    end

    if nargin < 7 || isempty(waves)
        waves = 0; % waves = 1 for wave propagators in stimulus and relay nuclei
    end
    
    if nargin < 6 || isempty(fs)
        fs = 10000;
    end
    
    if nargin < 5 || isempty(grid_edge)
        grid_edge = 12;
    else
        if mod(grid_edge,2)
            error('The number of grid edge nodes should be even to correctly align the FFT components')
        end
    end

    if nargin < 4 || isempty(int_time)
        int_time = 15; % Total integration time (s)
    end
    
    if nargin < 3 || isempty(firemode)
        firemode = [0 0 0 0];
    else
        firemode = [firemode(1) firemode(1) firemode(2:3)]; % Add an entry for inhibitory
    end
    
    if nargin < 2 || isempty(file_id)
        file_id = 1;
    end
    
    if size(p.nus,1) ~= 1 || size(p.phia,1) ~= 1
        error('Input point must have exactly one set of nuab and phia')
    end
	
    if p.Lx ~= p.Ly
        error('Rectangular grids are not currently supported using this interface. Please write your conf file manually');
    end

    if any([~isempty(p.spatial_gab) ~isempty(p.spatial_gabcd) ~isempty(p.spatial_alpha) ~isempty(p.spatial_beta) ~isempty(p.spatial_xyz)])
	   warning('Spatial variations are present that cannot be modeled in NeuroField. These will be ignored');
    end

    % Initialize numerical solver parameters
    deltat = 1/fs;
	deltax= p.Lx/grid_edge; % Assume square grid
	
	% CFL condition
	v = p.gammae*p.re; % Axonal velocity
    if deltat >= deltax/v % If the timestep is larger than the time required to propagate deltax
        fprintf(2,'Modifying deltat to satisfy CFL stability condition\n');
        deltat = deltax/v * 0.9;
        fprintf(1,'Deltat = %.10g\n',deltat);
    end
       
    confname = sprintf('%s_%d',fprefix,file_id);
    write_nf();

    if nargout > 0
        varargout{1} = nf.run(confname);
        varargout{1}.params = p.copy();
    else
        nf.run(confname);
    end
    
    if nargout > 3
        [f,P,V] = nf.spatial_spectrum(varargout{1},'propag.1.phi',[],8,1); % kmax=4, 8 segments, spatial filtering=1
        varargout{4} = V;
    elseif nargout > 1
        [f,P] = nf.spatial_spectrum(varargout{1},'propag.1.phi',[],8,1); % kmax=4, 8 segments, spatial filtering=1
    end
    
    if nargout > 1    
        varargout{2} = f;
        varargout{3} = P;
        figure
        loglog(f,P,'b');
        xlabel('Frequency (Hz)');
        ylabel('Power (arbitrary)');
        try
            title(sprintf('X: %.3f Y: %.3f Z: %.3f',p.xyz(1),p.xyz(2),p.xyz(3)));
        end
        set(gca,'XLim',[1 45]);

        [f_a,P_a] = p.spectrum;
        P_a = interp1(f_a,P_a,f,'pchip','extrap');
        hold on
        loglog(f,P_a,'r--');
        legend('NeuroField','Analytic');
        hold off

    end


    function write_nf()
        fid = fopen(sprintf('%s.conf',confname),'w');
        fprintf(fid,'EIRS model, automatically generated with nf.eirs.m\n');
        conn_matrix = [1 2 0 3 0; 4 5 0 6 0; 7 0 0 8 0; 9 0 10 0 11; 0 0 0 0 0].';
        labels = {'Excitatory','Inhibitory','Reticular','Relay','Stimulus'};
        phivals = [p.phia(1) p.phia(1) p.phia(2) p.phia(3)];
        n_dendrites = [3 3 2 3]; % Number of dendrites on each population
        id_map = [1 2 3 1 2 3 7 8 4 5 6]; % Map indices from 1-8 in a point struct to the 11 separate connections
        counter = 1; % A counter to keep track of the dendrites as they are added

        [a,b] = find(conn_matrix);
        for j = 1:length(a)
            fprintf(fid,'Connection %d - %s -> %s\n',j,labels{a(j)},labels{b(j)});
        end
        fprintf(fid,'\n');


        fprintf(fid,'Time: %.10g Deltat: %.10g\n',int_time,deltat);
        fprintf(fid,'Nodes: %i\n',grid_edge^2);
                         
        fprintf(fid,'\n');
                
        fprintf(fid,'    Connection matrix:\n');
        fprintf(fid,'From:  1  2  3  4  5\n');
        fprintf(fid,'To 1:  1  2  0  3  0\n');
        fprintf(fid,'To 2:  4  5  0  6  0\n');
        fprintf(fid,'To 3:  7  0  0  8  0\n');
        fprintf(fid,'To 4:  9  0  10 0  11\n');
        fprintf(fid,'To 5:  0  0  0  0  0 \n\n');


        for j = 1:4 % For each population
            fprintf(fid,'Population %d: %s\n',j,labels{j});
            fprintf(fid,'Length: %.10g\n',p.Lx);
            fprintf(fid,'Q: %.10g\n',phivals(j));
            if firemode(j) == 1 % If population has been linearized
                v0 = model.sinv(phivals(j),p);
                a = model.rho1(phivals(j),p);
                b = phivals(j) - a*v0;
                fprintf(fid,'Firing: Linear - a: %.10g b: %.10g\n',a,b);
                fprintf(1,'Linear %s\n',labels{j});
            elseif firemode(j) == 2 % If population has been quadraticized
                v0 = model.sinv(phivals(j),p);
                a = model.rho2(phivals(j),p)/2;
                b = model.rho1(phivals(j),p) - 2*a*v0;
                c = phivals(j) - model.rho1(phivals(j),p)*v0 + a*v0.^2;
                fprintf(fid,'Firing: Quadratic - a: %.10g b: %.10g c: %.10g\n',a,b,c);
                fprintf(1,'Quadratic %s\n',labels{j});
            elseif firemode(j) == 3 % Cubic population
                v0 = model.sinv(phivals(j),p);
                rho3 = @(Q,p) 6*Q^2*model.rho1(Q,p)/p.sigma^2/p.qmax^2 - 6*Q*model.rho1(Q,p)/p.sigma^2/p.qmax + model.rho1(Q,p)/p.sigma^2;
                r3 = rho3(phivals(j),p);
                r2 = model.rho2(phivals(j),p);
                r1 = model.rho1(phivals(j),p);
                d =  - (r3*v0^3)/6 + (r2*v0^2)/2 - r1*v0 + phivals(j);
                c =  (r3*v0^2)/2 - r2*v0 + r1;
                b = r2/2 - (r3*v0)/2;
                a = r3/6;
                fprintf(fid,'Firing: Cubic - a: %.10g b: %.10g c: %.10g d: %.10g\n',a,b,c,d);
                fprintf(1,'Cubic %s\n',labels{j});
            else
                fprintf(fid,'Firing: Sigmoid - Theta: %.10g Sigma: %.10g Qmax: %.10g\n',p.theta,p.sigma,p.qmax);
            end
            
            for k = 1:n_dendrites(j)
                fprintf(fid,' Dendrite %d: alpha: %.10g beta: %.10g\n',counter,p.alpha(id_map(counter)),p.beta(id_map(counter)));
                counter = counter + 1;
            end
            fprintf(fid,'\n');
        end

        fprintf(fid,'Population 5: Stimulation\n');
        fprintf(fid,'Length: %.10g\n',p.Lx);
        if isempty(ranseed)
            fprintf(fid,' Stimulus: White - Onset: 0 Mean: 1 Psd: %.10g\n',p.phin);
        else
            fprintf(fid,' Stimulus: White - Onset: 0 Ranseed: %d Mean: 1 Psd: %.10g\n',ranseed,p.phin);
        end
        x = deltax*(0:grid_edge-1);
        y = deltax*(0:grid_edge-1);
        [grid_x,grid_y] = meshgrid(x,y);

        if utils.isfunction(p.spatial_t0)
            % Vectors for the distance in each direction
            spec_t0 = p.spatial_t0(grid_x,grid_y);
            taustr = sprintf('%.10g ',spec_t0/2);
            fprintf('Spatial variations active\n');
        else
            taustr = sprintf('%.10g ',p.t0/2);
        end

        fprintf(fid,'\n');
        if waves
            error('placeholder')
        else
            fprintf(fid,'Propag 1: Wave - Tau: %.10g Range: %.10g gamma: %.10g\n',0,p.re,p.gammae);
            fprintf(fid,'Propag 2: Map - Tau: %.10g\n',0);   
            fprintf(fid,'Propag 3: Map - Tau: %s\n',taustr);   
            fprintf(fid,'Propag 4: Wave - Tau: %.10g Range: %.10g gamma: %.10g\n',0,p.re,p.gammae);
            fprintf(fid,'Propag 5:  Map - Tau: %.10g\n',0);  
            fprintf(fid,'Propag 6:  Map - Tau: %s\n',taustr); 
            fprintf(fid,'Propag 7: Wave - Tau: %s Range: %.10g gamma: %.10g\n',taustr,p.re,p.gammae);
            fprintf(fid,'Propag 8:  Map - Tau: %.10g\n',0);  
            fprintf(fid,'Propag 9: Wave - Tau: %s Range: %.10g gamma: %.10g\n',taustr,p.re,p.gammae);
            fprintf(fid,'Propag 10: Map - Tau: %.10g\n',0);  
            fprintf(fid,'Propag 11: Map - Tau: %.10g\n',0); % Stimulus wave
        end

        fprintf(fid,'\n');
        for j = 1:length(id_map)
            if ~isempty(p.spatial_nus) && utils.isfunction(p.spatial_nus{id_map(j)})
                couple_values = p.spatial_nus{id_map(j)}(grid_x,grid_y);
                fprintf('Spatial variations in Couple %d\n',j);
            else
                couple_values = p.nus(id_map(j));
            end
            couplestr = sprintf('%.10g ',couple_values);

            fprintf(fid,'Couple %d:  Map - nu: %s\n',j,couplestr);
        end
        fprintf(fid,'\n');

        fprintf(fid,'Output: Node: All Start: 5 Interval: 0.5e-2 \n');
        %fprintf(fid,'Output: Node: All Start: 0  \n');
        fprintf(fid,'Population: 1\n');
        fprintf(fid,'Dendrite:  \n');
        fprintf(fid,'Propag: 1\n');
        fprintf(fid,'Couple:  \n');

        fclose(fid);
        
        fprintf(1,'Integration time: %d s sampled at %d Hz\nGrid size %dx%d, outputting all nodes\nSimulating...',int_time,1/deltat,grid_edge,grid_edge);
    end
end
