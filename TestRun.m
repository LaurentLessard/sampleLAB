
%% Runs the function sampleLAB for different delta values
delta = 10:30 %linspace(10,30,20);
n = zeros(length(delta),1);
for d = 1:length(delta)
    coord = sampleLAB(delta(d));
    n(d) = length(coord(:,1));
end
%%
plot(delta,n)
title("Tradeoff between number of colors obtained vs \delta") 
xlabel('\delta');ylabel('Number of colors');
savefig('Tradeoff.fig');    
    