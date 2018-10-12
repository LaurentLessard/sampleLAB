%
rng(1);
max_muR = [0.3 0.9; 0.3 0.3; 0.7 0.9; 0.9 0.3]; % max x and y coordinate from all 4 boxes.
min_muR = [0.1 0.7; 0.1 0.1; 0.7 0.7; 0.7 0.1 ] ;% max x and y coordinate from all 4 boxes.
alphaR = [0.01 0.04]; %[min max]
betaR = [0.01 0.04];  %[min max]
window_size = 30;
noise_strength = 3;
delta = [2.2 2];
fp = plot_fig(max_muR, min_muR, alphaR, betaR, delta, window_size, noise_strength,'');

function Fp = plot_fig(max_muR, min_muR, alphaR, betaR, delta, window_size, noise_strength, cmap)
    % max_muR       : 4x2 matrix specifying max (x,y) value for mean in every quadrant.
    % min_muR       : 4x2 matrix specifying max (x,y) value for mean in for every quadrant.
    % alphaR        : 1x2 vector with range of alpha [min max]
    % betaR         : 1x2 vector with range of alpha [min max]
    % delta         : 1x2 vector with values delta(1) > delta(2)
    % window_size   : scalar
    % noise_strength: scalar
    % cmap          : colormap
    
    x1 = 0:0.01:1;
    x2 = x1;
    n = 25;
    fig = tight_subplot(5,5,[.03 .02],[.03 .01],[.01 .01]);
    for i = 1:n
        [mu1,mu2,mu3,mu4] = get_means(max_muR, min_muR)
        [Sigma1, Sigma2, Sigma3, Sigma4] = get_covariance(alphaR,betaR)
        [f1,f2,f3,f4] = get_mvnpdf(mu1,mu2,mu3,mu4,Sigma1,Sigma2,Sigma3,Sigma4);
        [amp1, amp2, amp3, amp4] = get_amplitudes(delta);
        F = (amp1*f1) + (amp2*f2) + (amp3*f3) + (amp4*f4);
        % Generate Gaussian noise
        noise = (randn(size(F)));
        x = linspace(-2,2,window_size);
        ex = exp(-x.^2);
        H = ex'*ex;
        H = H / sum(sum(H));
        % Filter noise
        Y = filter2(H,noise,'same');
        Fp = reshape(F+(noise_strength*Y),size(F)); % Noise Strength
        %figure(8); surf(x1,x2,Fp);
        %subplot(5,5,i);%set(gcf, 'Position', get(0, 'Screensize'));
        axes(fig(i))
        contourf(x1,x2,Fp,8,'LineStyle','none');
        %daspect([1 1 1]);
        colormap(cmap);
%         frames(i) = getframe(gcf);
%         pause(0.5);
    end
%     video = VideoWriter('Plots.avi','Uncompressed AVI');
%     video.FrameRate = 2;
%     open(video);
%     writeVideo(video,frames);
%     close(video);
end

function [m1, m2, m3, m4] = get_means(max_muR, min_muR)
    for i = 1:4
        m = min_muR + (max_muR-min_muR).*(rand(size(max_muR)));
        A = m;
    end
    m1 = A(1,:);
    m2 = A(2,:);
    m3 = A(3,:);
    m4 = A(4,:);
end

function [sigma1, sigma2, sigma3, sigma4,theta] = get_covariance(alphaR,betaR)
    for i = 1:4 
        a = alphaR(1) + (alphaR(2)-alphaR(1))*rand;
        b = betaR(1) + (betaR(2)-betaR(1))*rand;
        theta(i) = 360.*rand;
        u = [cosd(theta(i)) sind(theta(i)); -sind(theta(i)) cosd(theta(i))];
        covar = u'*[a 0; 0 b]*u;
        if i>1
            A(:,:,i) = covar;
        else
            A = covar;
        end
    end
    sigma1 = A(:,:,1);
    sigma2 = A(:,:,2);
    sigma3 = A(:,:,3);
    sigma4 = A(:,:,4);
end

function [f1,f2,f3,f4] = get_mvnpdf(mu1,mu2,mu3,mu4,Sigma1,Sigma2,Sigma3,Sigma4)
    x1 = 0:0.01:1;
    x2 = x1;
    [X1,X2] = meshgrid(x1,x2);
    f1 = mvnpdf([X1(:) X2(:)],mu1,Sigma1);
    f1 = reshape(f1,length(x2),length(x1));
    f2 = mvnpdf([X1(:) X2(:)],mu2,Sigma2);
    f2 = reshape(f2,length(x2),length(x1));
    f3 = mvnpdf([X1(:) X2(:)],mu3,Sigma3);
    f3 = reshape(f3,length(x2),length(x1));
    f4 = mvnpdf([X1(:) X2(:)],mu4,Sigma4);
    f4 = reshape(f4,length(x2),length(x1));
end
function [amp1, amp2, amp3, amp4] = get_amplitudes(delta)
    amp1 = (0.5 + (0.2*rand))*delta(1); % amp1 lies within 50-70% of delta1
    amp2 = delta(1) - amp1;
    % amp3 is some percent of delta2
    while 1
        amp3 = rand*delta(2);
        amp4 = delta(2) - amp3;
        % Check for the inequalities
        if(amp1 + amp4 > amp3 + amp2)
            break
        end
    end
end