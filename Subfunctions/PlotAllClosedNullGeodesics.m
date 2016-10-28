% function PlotAllClosedNullGeodesics(x1Psol,x2Psol,x1_g,x2_g,lamV,lam2)

% Input arguments:
    % lamV : Desired set of \lambda values
    % phi0 : initial \phi value (cf. Fig. 2 of [1])
    % CGij : ij entries of the CG strain tensor 
    % x1_g  : x1 component of the spatial grid 
    % x2_g  : x2 component of the spatial grid 

%--------------------------------------------------------------------------
% Author: Mattia Serra  serram@ethz.ch
% http://www.zfm.ethz.ch/~serra/
%--------------------------------------------------------------------------
function PlotAllClosedNullGeodesics(x1Psol,x2Psol,x1_g,x2_g,lamV,lam2)
    

    % Colormap encoding different \lambda values 
    cmap = jet(length(lamV));
    % Plot properties 
    AxthicksFnt = 15;
    fontsizeaxlab = 15;
    
    if ~isempty(x1Psol) % If there are closed null-geodesics 

     % Initialize the figure with the FTLE plot 
    figure('units','normalized','outerposition',[0 0 .5 .5]);
    imagesc(x1_g,x2_g,log(lam2)/30/2);shading interp
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


    for kkmuv=1:1:length(lamV)
    xxapp=x1Psol{kkmuv};
    yyapp=x2Psol{kkmuv};
    if ~isempty(xxapp)
        for kkc=1:size(xxapp,2)
            xlc=xxapp(~isnan(xxapp(:,kkc)),kkc);
            ylc=yyapp(~isnan(yyapp(:,kkc)),kkc);         
            hold on
            plot(xlc,ylc,'color',cmap(kkmuv,:),'linewidth',2.5)             
        end
    end
    end
    axis equal tight

    %add a second colorbar for the \lambda values
    ax1=gca;
    ax1_pos = ax1.Position; 
    ax2 = axes('Position',ax1_pos,...
               'XAxisLocation','bottom',...
               'YAxisLocation','left',...
               'Color','none');
    hhF2 = colorbar(ax2,'eastOutside')
    hhF2.FontSize=AxthicksFnt;
    set(get(hhF2,'xlabel'),'string','$$\lambda$$','Interpreter','latex','FontWeight','normal');
    colormap(ax2,'jet')
    hhF2.Ticks=linspace(0,1,3);
    hhF2.XTickLabel={'0.9';'1';'1.1'};
    set(ax2,'xtick',[])
    set(ax2,'ytick',[])
    set(ax2, 'visible', 'off') ;
    end
end