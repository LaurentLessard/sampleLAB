
%% Runs the function sampleLAB for different delta values
delta = linspace(10,40,10);
n = zeros(length(delta),1);
for d = 1:length(delta)
    coord = sampleLAB(delta(d));
    n(d) = length(coord(:,1))
end
%%
plot(delta,n)
title("Tradeoff between number of colors obtained vs delta")  
%%
    

    