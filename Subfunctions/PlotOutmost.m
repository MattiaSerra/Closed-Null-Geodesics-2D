% function PlotOutmost(xLcOutM,yLcOutM,LamLcOutM,lamV,x_g,y_g,lam2)

% Input arguments:
    % xLcOutM                    : x-component of the outermost closed null-geodesics 
    % yLcOutM                    : x-component of the outermost closed null-geodesics 
    % LamLcOutM                  : \lambda values of the outermost closed null-geodesics 
    % lamV,lamV,x_g,y_g,lam2     : see step 1

%--------------------------------------------------------------------------
% Author: Mattia Serra  serram@ethz.ch
% http://www.zfm.ethz.ch/~serra/
%--------------------------------------------------------------------------
function PlotOutmost(xLcOutM,yLcOutM,LamLcOutM,lamV,x_g,y_g,lam2)

     if ~isempty(LamLcOutM) % If there are closed null-geodesics 
        % Colormap encoding different \lambda values 
        cmap = jet(length(lamV));

        % Plot properties 
        AxthicksFnt = 15;
        fontsizeaxlab = 15;

        % Initialize the figure with the FTLE plot 
        figure('units','normalized','outerposition',[0 0 .5 .5]);
        imagesc(x_g,y_g,log(lam2)/30/2);shading interp
        set(gca,'FontSize',AxthicksFnt,'fontWeight','normal')
        hold on
        set(gca,'YDir','normal')
        set(gcf,'color','w');
        axis equal
        xlabel('$$Lon [^{\circ}]$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
        ylabel('$$Lat [^{\circ}]$$','Interpreter','latex','FontWeight','bold','FontSize',fontsizeaxlab);
        axis equal tight 
        colormap(gca,'gray')
        hhF=colorbar(gca);
        hhF.Location='westOutside';
        hhF.FontSize=fontsizeaxlab;
        set(get(hhF,'xlabel'),'string','$$FTLE$$','Interpreter','latex','FontWeight','normal');


        % Plot outermost Closed null-geodesics 
        for kkmuv=1:1:length(LamLcOutM)
                Lamidx=find(lamV==LamLcOutM(kkmuv));
                xlc=xLcOutM{kkmuv};
                ylc=yLcOutM{kkmuv};
                hold on
                plot(xlc,ylc,'color',cmap(Lamidx,:),'linewidth',2.5);
        end
        axis equal tight

        % Add a second colorbar encoding the different \lambda values 

        ax1 = gca;
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
                   'XAxisLocation','bottom',...
                   'YAxisLocation','left',...
                   'Color','none');
        hhF2 = colorbar(ax2,'eastOutside');
        hhF2.FontSize = AxthicksFnt;
        set(get(hhF2,'xlabel'),'string','$$\lambda$$','Interpreter','latex','FontWeight','normal');
        colormap(ax2,'jet')
        hhF2.Ticks=linspace(0,1,3);
        hhF2.XTickLabel={'0.9';'1';'1.1'};
        set(ax2,'xtick',[])
        set(ax2,'ytick',[])
        set(ax2, 'visible', 'off') ;
     end
     
end