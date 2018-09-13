function coord = sampleLAB(delta)
    %delta = 10
    % Generate delta spaced points for L (0-100),a(-100 to 100) and b (-100 to 100)
    L = 0:delta:100 % length = N
    a = -100:delta:100 % length = 2N-1
    b = -100:delta:100 % length = 2N-1
    lab = (combvec(L,a,b))'
    % Every row of the matrix 'lab' is a point in the color space
    % Validate coordinates to exist in RGB space
    rgb = lab2rgb(lab)
    [row,col] = find(rgb > 1 | rgb < 0)
    indices = unique(row)
    rgb(indices,:) = [] 
    % Get rid of coordinates which do not exist in RGB space
    lab(indices,:) = []
    out = colorconvert(lab,'Lab','D65')
    xyY = [out.x out.y out.Y]
    coord = xyY
    scatter3(xyY(:,1),xyY(:,2),xyY(:,3),50,rgb,'filled')
end
