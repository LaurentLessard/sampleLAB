function lab = sampleLAB(delta)
    % sampleLAB(delta): Generates a grid of points uniformly spaced in LAB space 
    % which can be mapped to RGB space.
    % Input : delta -- TODO
    % Output : Nx3 matrix of coordinates in xyY space.
    close all; clc;
    % Generate delta spaced points for L (0-100),a(-100 to 100) and b (-100 to 100)
    L = 0:delta:100 ;% length = N
    a = [0-delta:-delta:-100 0:delta:100];
    b = a;
    lab = (combvec(L,a,b))';
    %Plot
    figure;  subplot(1,2,1);set(gcf, 'Position', get(0, 'Screensize'));
    scatter3(lab(:,2),lab(:,3),lab(:,1),10,'k');
    hold on;
    % Every row of the matrix 'lab' is a point in the color space
    % Validate coordinates to exist in RGB space
    rgb = lab2rgb(lab);
    [row,~] = find(rgb > 1 | rgb < 0);
    indices = unique(row);
    % Get rid of coordinates which do not exist in RGB space
    rgb(indices,:) = [] ;
    lab(indices,:) = [];
    scatter3(lab(:,2),lab(:,3),lab(:,1),80,rgb,'filled');
    daspect([1 1 1]); axis([-100 100 -100 100 0 100]); grid on; 
    title(['LAB space with uniformly distributed points for \delta = ',num2str(delta)]);
    legend(["Uniformly distributed points","Selected points existing in RBG space"],'Location','southeast');
    xlabel('a');ylabel('b'); zlabel('L');
    hold off;
    out = colorconvert(lab,'Lab','D65');
    xyY = [out.x out.y out.Y];
    lab = xyY;
    subplot(1,2,2);set(gcf, 'Position', get(0, 'Screensize'));
    scatter3(xyY(:,1),xyY(:,2),xyY(:,3),70,rgb,'filled');
    axis([0 1 0 1 0 100]); daspect([0.1 0.1 10]);
    title(['Selected colors in xyY space for \delta = ',num2str(delta)]);
    xlabel('x');ylabel('y'); zlabel('Y');
    savefig(strcat('fig_delta_',num2str(delta),'.fig'));
end
