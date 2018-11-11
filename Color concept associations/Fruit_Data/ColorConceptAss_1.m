% For 1 concept - Later generalize
% Using first 15 observation data
obs = 20 %change to 54 % given 54 data files
A = zeros(obs,58);
B = zeros(obs,1);

% get data for a particular concept
for i = 1:obs
    %A(:,i) = xlsread(strcat(num2str(i),"-ColorConceptAssocFr.xls"),'D2:D59') % for concept mango)
    T = table2array(readtable(strcat(num2str(i),"-ColorConceptAssocFr.xls"),'Range','D2:D59','ReadVariableNames',false));
    % A: (#obs x #colors) Matrix of each observant's color weights. # Colors = 58 
    A(i,:) = T';
    [val, ind] = max(T);
    % B vector of the colors each observant thought best suited the concept
    B(i,1) = ind; 
end
p = mean(A);
[v,ind] = max(p);
ind
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the corelation coefficient for every color given a concept
figure;
scatter(1:size(A,1),A(:,58));
