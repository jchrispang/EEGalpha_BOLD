function k_explorer(f,P,Kx,Ky)
    if nargin < 3
        [Kx,Ky] = meshgrid(1:size(P,1),1:size(P,2));
    end
    
    handles.figure1 = gcf;
    handles.f = f;
    handles.P = P;
    handles.Kx = Kx;
    handles.Ky = Ky;
    
    subplot(1,2,1);
    axis equal
    handles.s = patch([0],[0],'r');
    set(handles.s,'FaceColor','r');
    hold on
    handles.l = plot([0 0 NaN 0 0],[0 0 NaN 0 0],'r--');
    hold off
    box on
    
    set(gca,'XLim',[1,size(Kx,2)],'YLim',[1,size(Ky,1)],'XTick',[1:size(Kx,2)],'YTick',[1:size(Ky,1)],'XTickLabel',Kx(1,:),'YTickLabel',Ky(:,1))
    grid on
    handles.a1 = gca;

    set(handles.a1,'ButtonDownFcn',@(a,b,c) refreshplot(guidata(a)));
    
    subplot(1,2,2);
    loglog(2,2);
    
    guidata(handles.figure1,handles);
    
function refreshplot(handles)
    pos = get(handles.a1,'CurrentPoint');
    xy = round(pos(1,1:2));
    set(handles.s,'XData',[xy(1)-0.5 xy(1)-0.5 xy(1)+0.5 xy(1)+0.5],'YData',[xy(2)+0.5,xy(2)-0.5,xy(2)-0.5,xy(2)+0.5]);
    set(handles.l,'XData',[1 xy(1) NaN xy(1) xy(1)],'YData',[xy(2) xy(2) NaN 1 xy(2)]);
    subplot(1,2,2);
    loglog(handles.f,squeeze(handles.P(xy(1),xy(2),:)));
    title(sprintf('K_x = %.2f, K_y = %.2f',handles.Kx(1,xy(1)) ,handles.Ky(xy(2),1)));
    set(gca,'XLim',[1 45]);

