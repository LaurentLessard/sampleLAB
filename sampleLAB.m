function [lab, rgb, xyY] = sampleLAB(delta, rot)
    % sampleLAB(delta,rot): Generates a grid of points uniformly spaced in LAB space 
    % which can be mapped to RGB space.
    % Input : delta (uniform distance between points) , rot (rotation angle)
    % Output : Nx3 matrix of coordinates in xyY space.
    
    clc;
    % Generate delta spaced points for L (0-100),a(-100 to 100) and b (-100 to 100)
    L = 0:delta:100 ;% length = N
    a = [0-delta:-delta:-100*sqrt(2) 0:delta:100*sqrt(2)];
    b = a;
    
    lab = (combvec(L,a,b))';
    labOrig = lab; 
    lab = lab * [ 1 0 0; 0 cosd(rot) sind(rot); 0 -sind(rot) cosd(rot)];
    %Plot
    figure;  
    subplot(1,2,1);set(gcf, 'Position', get(0, 'Screensize'));
    scatter3(lab(:,2),lab(:,3),lab(:,1),10,'k');
    hold on;
    
    % Every row of the matrix 'lab' is a point in the color space
    % Validate coordinates to exist in RGB space
    rgb = lab2rgb(lab);
    [row,~] = find(rgb > 1+eps | rgb < 0-eps);
    indices = unique(row);
    
    % Get rid of coordinates which do not exist in RGB space
    rgb(indices,:) = [] ;
    lab(indices,:) = [];
    labOrig(indices,:) = [];
    
    % Plot the selected points
    scatter3(lab(:,2),lab(:,3),lab(:,1),80,rgb,'filled','MarkerEdgeColor',[0 0 0]);
    daspect([1 1 1]); axis([-100 100 -100 100 0 100]); grid on;
    title(['LAB space with uniformly distributed points for \delta = ',num2str(delta),' and rotation angle =',num2str(rot)]);
    legend(["Uniformly distributed points","Selected points existing in RBG space"],'Location','southeast');
    xlabel('a');ylabel('b'); zlabel('L');
    hold off;
    
    % Convert Lab to xyY
    out = colorconvert(lab,'Lab','D65');
    xyY = [out.x out.y out.Y];
    
    % Plot colors in xyY space
    subplot(1,2,2);set(gcf, 'Position', get(0, 'Screensize'));
    scatter3(xyY(:,1),xyY(:,2),xyY(:,3),70,rgb,'filled','MarkerEdgeColor',[0 0 0]);
    axis([0 1 0 1 0 100]); daspect([0.1 0.1 10]);
    title(['Selected colors in xyY space for \delta = ',num2str(delta),' and rotation angle =',num2str(rot)]);
    xlabel('x');ylabel('y'); zlabel('Y');
    txt = ['# Colors: ',num2str(length(rgb(:,1)))];
    col = length(rgb(:,1));
    text(0.2,1,100,txt);
    % Save figure
    savefig(strcat('output2/fig_delta_',num2str(delta),'.fig'));
    csvwrite('Lab.csv',lab)
    csvwrite('RGB.csv',rgb)
    csvwrite('xyY.csv',xyY)
    csvwrite('LabOrig.csv',labOrig)
end
