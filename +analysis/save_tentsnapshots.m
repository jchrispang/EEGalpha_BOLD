function save_tentsnapshots(fname)
    % This function saves 4 views of the output from bin_3d
  	%xray = get(gcf,'UserData');
	origview = get(gca,'View');

	set(gca,'View',[117 12]) % iso
	utils.pfig(sprintf('%s_iso',fname),'png')

	set(gca,'View',[90 0]) % front
	%xray(1)
	utils.pfig(sprintf('%s_front',fname),'png')
	%xray(0)

	set(gca,'View',[180 0])
	%xray(1)
	utils.pfig(sprintf('%s_side',fname),'png')
	%xray(0)

	set(gca,'View',[0 90]) % top
	%xray(1)
	utils.pfig(sprintf('%s_top',fname),'png')
	%xray(0)

	set(gca,'View',origview)

	saveas(gcf,sprintf('~/Desktop/%s.fig',fname),'fig')
