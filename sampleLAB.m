%%
function lab = sampleLAB(delta)
    % sampleLAB(delta): Generates a grid of points uniformly spaced in LAB space 
    % which can be mapped to RGB space.
    % Input : delta --TODO
    % Output : Nx3 matrix of coordinates in RGB space.
    close all; clc;
    % Generate delta spaced points for L (0-100),a(-100 to 100) and b (-100 to 100)
    L = 0:delta:100 ;% length = N
    a = -100:delta:100; % length = 2N-1
    b = -100:delta:100; % length = 2N-1
    lab = (combvec(L,a,b))';
    figure; set(gcf, 'Position', get(0, 'Screensize'));
    scatter3(lab(:,1),lab(:,2),lab(:,3),30,'green');
    hold on;
    % Every row of the matrix 'lab' is a point in the color space
    % Validate coordinates to exist in RGB space
    rgb = lab2rgb(lab);
    [row,col] = find(rgb > 1 | rgb < 0);
    indices = unique(row);
    % Get rid of coordinates which do not exist in RGB space
    rgb(indices,:) = [] ;
    lab(indices,:) = [];
    scatter3(lab(:,1),lab(:,2),lab(:,3),50,'r','filled');
    title("LAB space with uniformly distributeed coordinates. ");
    legend(["Uniformly distributed coordinates","Selected points existing in RBG space"]);
    savefig('Fig1.fig');
    out = colorconvert(lab,'Lab','D65');
    xyY = [out.x out.y out.Y];
    coord = xyY;
    figure; set(gcf, 'Position', get(0, 'Screensize'));
    scatter3(xyY(:,1),xyY(:,2),xyY(:,3),50,rgb,'filled');
    savefig('Fig2.fig'); title("xyY space");
end
