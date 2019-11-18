function flipbook(prefix)
	page_size = [15 8]; % Number of plots per page
	paper_size = [29 21]*2;
	xlimits = [1 45];
	ylimits = [1e-4 2];

 
    [main_figure,h] = initialize_figure(page_size,paper_size);

    
    files = dir(sprintf('sample_spectra/%s',prefix));
    files = files(4:end);
    count = 1;
    for j = 0:length(files)-1
    	    fname = files(j+1).name;
		    S = load(sprintf('sample_spectra/%s/%s',prefix,fname),'f','P','scaling');

			f = S.f;
			P = S.P;

			if mod(j,prod(page_size))==0 && j~=0			    
			    set(h(1),'YTick',ylimits(1),'YTickLabel',sprintf('%.0e',ylimits(1)));
			    set(h(1+page_size(2)),'YTick',ylimits(2),'YTickLabel',sprintf('%.0e',ylimits(2)));          
			    for col = 1:page_size(2)
                    index = col;
                    set(h(col),'XTick',xlimits,'XTickLabel',xlimits);
			    end 
			    
			    
			    pfig(sprintf('book_%s_p%i',prefix,count),'eps',main_figure);

			    fprintf(1,'Printed %i of %i\n',j+1,length(files)+1);
			    count = count + 1;

			    [main_figure,h] = initialize_figure(page_size,paper_size);

			end
			
			this_axis = h(1+mod(j,prod(page_size)));
			axes(this_axis)
			loglog(this_axis,f,P,'LineWidth',2);
			hold(this_axis,'on');
			plot(this_axis,[10 10],ylimits,'r--'); % 10 Hz line
			plot(this_axis,[20 20],ylimits,'r--'); % 20 Hz line
			
			set(this_axis,'XLim',xlimits,'YLim',ylimits)
            set(this_axis,'XTickLabel',[],'YTickLabel',[]);
            text(0.02,0.1,strrep(fname,'_','\_'),'HorizontalAlignment','left','VerticalAlignment','mid','FontSize',8,'Units','normalized')
            
	end
	pfig(sprintf('book_%s_p%i',prefix,count));
end

function [main_figure,h] = initialize_figure(page_size,paper_size)
    close all
	main_figure = figure;
	
	left = linspace(1.5,paper_size(2),page_size(2)+1);
    bottom = linspace(1,paper_size(1),page_size(1)+1);
	width = paper_size(2)/page_size(2);
	height = paper_size(1)/page_size(1);
	   
    set(main_figure,'PaperPositionMode','auto')
    set(main_figure,'Units','centimeters');
    set(main_figure,'Position',[0   0 paper_size(2)   paper_size(1)])

    for row = 1:page_size(1)
        for col = 1:page_size(2)
            index = (row-1)*page_size(2)+col;
            h(index) = axes('Units','centimeters','position',[left(col) bottom(row)  width*0.95  height*0.95]); 
        end
    end
end			    
	
