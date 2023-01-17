function plot_boundaries_poles(IM,Xp,Yp,phirot,IPoles,path_data)

nSlices = size(IM,4);

figure('Visible','off')

for i=1:nSlices
    I = max(IM(:,:,:,i),[],3);
%     I  = imrotate(I,phirot);

    xp = round(Xp(:,i));
    yp = round(Yp(:,i));
%     [xp,yp] = rotatePoints(phirot-pi/2,xp,yp);
    if nSlices > 1
        subaxis(3,4,i)
    end
    imshow(I,[]); hold on;
    plot(xp,yp,'.','Color',[0 0.4470 0.7410])
    plot(xp(IPoles(i,1)),yp(IPoles(i,1)),'*','MarkerSize',12,'color','red')
    plot(xp(IPoles(i,2)),yp(IPoles(i,2)),'*','MarkerSize',12,'color','yellow')
    text(xp(IPoles(1)),yp(IPoles(1)),'A','color','red','linewidth',10,'fontsize',12,'HorizontalAlignment','left') 

end


end