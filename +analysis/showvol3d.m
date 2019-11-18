function showvol3d(colour,transparency)
    % This function takes colour and transparency matrices and shows the plot
	figure(1)
	set(gcf,'Renderer','OpenGL')
	p = vol3d_alpha('Cdata',colour,'texture','3D','alphadata',log(transparency));
	vol3d_alpha(p);
	colorbar;
	p = vol3d_alpha('Cdata',colour,'texture','3D','alphadata',log(transparency));
	grid on
	
	
	p = vol3d_alpha('Cdata',colour,'texture','3D','alphadata',log(transparency));
	
	xyzlim = [0 -1 0 ;...
	          1  1 1.3];
	
	%caxis(min(min(min(colour>>0))),max(max(max(colour))));
	
			nticks = [5,11,5];
			xticks = linspace(xyzlim(1,1),xyzlim(2,1),nticks(1));
			yticks = linspace(xyzlim(1,2),xyzlim(2,2),nticks(2));
			zticks = linspace(xyzlim(1,3),xyzlim(2,3),nticks(3));

			xBinIndex = get(gca,'XLim');
			yBinIndex = get(gca,'YLim');
			zBinIndex = get(gca,'ZLim');
					
			set(gca,'XLim',[min(xBinIndex),max(xBinIndex)]);
			set(gca,'YLim',[min(yBinIndex),max(yBinIndex)]);
			set(gca,'ZLim',[min(zBinIndex),max(zBinIndex)]);
			set(gca,'XTick',linspace(min(xBinIndex),max(xBinIndex),nticks(1)));
			set(gca,'YTick',linspace(min(yBinIndex),max(yBinIndex),nticks(2)));
			set(gca,'ZTick',linspace(min(zBinIndex),max(zBinIndex),nticks(3)));
			set(gca,'XTickLabel',xticks,'YTickLabel',yticks,'ZTickLabel',zticks);
			
			xlabel('X')
			ylabel('Y')
			zlabel('Z')
