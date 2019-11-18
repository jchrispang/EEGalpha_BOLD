function [t,V] = erp(p)
	% Return the erp
	% Adapted from Rhys' forwardmodelORIGINAL, originally
	% written by Cliff

	maxfrequency = 50; % Choose maximum frequency in Hz
	maxtime      = 1; % Maximum time in seconds
	lnorm        = 1;
	tOffset      = 0.05;
	twImp        = 0.01;
	rOffset      = 0.1;
	rwImp        = 0.01;
	v            = 10.0;

	% Define other parameters used in the integration etc.
	delta_t         = 0.001;
	NPoints         = 10000;
	antiwrap        = 1; % Note: in erpfit this is one
	NPoints         = NPoints*antiwrap;
	eye             = complex(0,1);
	oversamp        = 1;
	maxFactorQr     = 2.0;
	maxFactorGauss  = 3.0;
	maxFactorBessel = 12.0;
	kreSteps        = 20;
	Tres            = delta_t/oversamp;
	timesmall       = (0:NPoints)*Tres*oversamp;
	omega           = zeros(NPoints,1);

	for aa=NPoints-1:-1:0
	    omega(aa+1)=2*pi*aa*(1.0/NPoints/Tres);
	end
	erpsmall=zeros(fix(NPoints/oversamp),1);


	% Decide where to place the cutoff -- only used for sum-over-modes
	% approach, not integral approach
	cutoff(1)=abs(maxFactorQr*4.0);
	cutoff(2)=abs(maxFactorGauss*(p.re/rwImp));
	cutoff(3)=abs(maxFactorBessel*(p.re/rOffset));
	kreMax=min(cutoff);

	j0vals=zeros(kreSteps,1);
	for n=0:kreSteps-1;
	    j0vals(n+1)=besselj(0,n*(kreMax/kreSteps)*rOffset/p.re);
	end

	Lee=1.0./((1-eye.*omega./p.alpha(1)).*(1-eye.*omega./p.beta(2)));
	Lii=Lee;
	Lre=Lee;
	Lsr=Lee;
	Les=Lee;
	Lsn=Lee;
	Lrs=Lee;
	lag1 = cos(omega*p.t0/2)+eye*sin(omega*p.t0/2);
	lag2 = cos(omega*p.t0)+eye*sin(omega*p.t0);
	lag3 = cos(omega*tOffset)+eye*sin(omega*tOffset);
	lowpass_c = Lee./(1-p.gabcd(2)*Lii);
	lowpass_t = 1./(1-p.gabcd(5).*Lsr.*Lrs);
	gain1 = Lsn.*lowpass_t.*lag1;
	gain2 = p.gabcd(3)+p.gabcd(4).*Lre;
	gain3 = p.gabcd(1)+Les.*lowpass_t.*gain2.*lag2;
	q2r2 = (1-eye.*omega./p.gammae).*(1-eye.*omega./p.gammae)-lowpass_c.*gain3;
	gauss1D = exp(-omega.*omega.*twImp.*twImp./2);

	k_int=0.0;
	for n=0:kreSteps-1
	    kr=n*kreMax/kreSteps;
	    k2r2=kr*kr;
	    k_int1=kr./(k2r2+q2r2);
	    k_int2=exp(-(kr*rwImp/(2*p.re))^2);%-exp(-k2r2*(2*rwImp/re/2)^2) % Change the stimulus!
	    % k_int2=k2r2*exp(-(kr*rwImp/(2*re))^2) % Mexican hat stimulus
	    k_int3=j0vals(n+1);
	    k_int=k_int+k_int2.*k_int3.*k_int1;
	end
	k_int=kreMax/kreSteps/2/pi/p.re/p.re.*k_int;

	tmp = k_int;
	tmp = gauss1D.*lag3.*lowpass_c.*gain1.*tmp;
	tmp = lnorm*tmp;

	Y=tmp;
	for ab=0:NPoints/2
	    Y(NPoints-ab)=conj(Y(ab+2)); % Fixed--should now exactly reproduce erpfit's results
	end

	y=fft(Y);
	y=y./Tres./NPoints;
	erpbig=real(y);
	for ab=1:NPoints/oversamp
	    cc=fix(timesmall(ab)/Tres+1);
	    erpsmall(ab)=erpbig(cc);
	end

	%% Save output data, and trim to only use maxtime and maxfreq points

	% ERP time series
	timeseries.x=(0:length(erpsmall)-1)'*Tres*oversamp; % Time
	timeseries.y=erpsmall; % Voltage as a function of time
	timepointstouse=(timeseries.x<=maxtime);
	timeseries.x=timeseries.x(timepointstouse);
	timeseries.y=timeseries.y(timepointstouse);

	% Write outputs
	t = timeseries.x;
	V = timeseries.y;