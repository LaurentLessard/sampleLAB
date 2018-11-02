% Required input for the function plot_fig
rng(1);
max_muR = [0.3 0.95; 0.2 0.3; 0.95 0.85; 0.95 0.3]; % max x and y coordinate from all 4 boxes.
min_muR = [0.2 0.9; 0.1 0.2; 0.85 0.8; 0.9 0.2 ] ;% max x and y coordinate from all 4 boxes.
alphaR = [0.05 0.07]; %[min max]
betaR = [0.05 0.07];  %[min max]
window_size = 25 ;%25 %10% 21%30 %28 % 25 %28;
noise_strength = 8 ;%10%2.2;
delta = [5 1];

figure(1)
clf;

fp = plot_fig(max_muR, min_muR, alphaR, betaR, delta, window_size, noise_strength);

function Fp = plot_fig(max_muR, min_muR, alphaR, betaR, delta, window_size, noise_strength)
    % max_muR       : 4x2 matrix specifying max (x,y) value for mean in every quadrant.
    % min_muR       : 4x2 matrix specifying max (x,y) value for mean in for every quadrant.
    % alphaR        : 1x2 vector with range of alpha [min max]
    % betaR         : 1x2 vector with range of alpha [min max]
    % delta         : 1x2 vector with values delta(1) > delta(2)
    % window_size   : scalar
    % noise_strength: scalar
    
    r = 1; % As flag for figure number in the row
    n = 40 ;  % n = Number of Replications * 2
    version = 0;
    x1 = linspace(0,1,200);
    x2 = linspace(0,1,100);
    %close all;
    %mkdir('outputImages');
    for i = 1:n
        if rem(i,2) == 1 % Left Heavy 
            [mu1,mu2,mu3,mu4] = get_means(max_muR, min_muR);
            [Sigma1, Sigma2, Sigma3, Sigma4] = get_covariance(alphaR,betaR);
            [f1,f2,f3,f4] = get_mvnpdf(mu1,mu2,mu3,mu4,Sigma1,Sigma2,Sigma3,Sigma4,x1,x2);
            [amp1, amp2, amp3, amp4] = get_amplitudes(delta, 'L');
            close all;
            %figure(fig_num); fig_num = fig_num + 1;
            version = version + 1;
            %mkdir('output',num2str(fold_num));
            param1 = 'r'; % for figure names
            param2 = 'l';
            r = 1;
        else  % Right heavy figures
            [amp1, amp2, amp3, amp4] = get_amplitudes(delta, 'R');
            param1 ='l';
            param2 = 'r';
        end
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
        
        % resize the image by downsampling (PICK AN EVEN NUMBER)
        N2 = 24;
        N1 = 12;
        Fp_sq = imresize(Fp,[N1 N2]);
        
        % randomize each half of the image
        q = 8;
        M = N1*N2/q;
        rix = randperm(M);
        for j = 1:q-1
            rix = [rix j*M+randperm(M)]; %#ok<AGROW>
        end
        Fp_sq = reshape(Fp_sq(rix),N1,N2);
        
        for c = 1:2
            if c == 1
                cmap = 'hot';col='hot';
            else
                cmap = viridis;col='vir';
            end
            % Each Fp to be in 2 color scales, 2 spatial arrangements
                h1 = figure ; %subplot(4,4,r); %figure;%; 
                r= r+1;
                contourf(Fp,8,'LineStyle','none'); colormap(h1,cmap);
                axis image
                axis off;%set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf,strcat('outputImages/',col,'_',param1,'_conc_light_',num2str(version),'.bmp'));
                
                h2 = figure ;% subplot(4,4,r); %figure;%subplot(4,4,r);
                r= r+1;
                imagesc(Fp_sq); colormap(h2,cmap);
                %set(h2,'YDir','normal')
                set(gca,'YDir','normal')
                axis image
                axis off;%set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf,strcat('outputImages/',col,'_',param1,'_grid_light_',num2str(version),'.bmp'));

                h3 = figure ;% subplot(4,4,r); %figure;%(4,4,r); 
                r= r+1;
                contourf(-Fp,8,'LineStyle','none');colormap(h3,cmap);
                axis image
                axis off;%set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf,strcat('outputImages/',col,'_',param2,'_conc_dark_',num2str(version),'.bmp'));

                h4 = figure ; %subplot(4,4,r); %figure;%subplot(4,4,r); 
                r= r+1;
                imagesc(-Fp_sq); colormap(h4,cmap);
                %set(h4,'YDir','normal')
                set(gca,'YDir','normal')
                axis image
                axis off;%set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf,strcat('outputImages/',col,'_',param2,'_grid_dark_',num2str(version),'.bmp'));
                
        end  
%        close all;
%         if rem(i,2) == 0
%             savefig(strcat('output/plot_',num2str(i/2),'.fig'));
%         end
    end
    
end

function [m1, m2, m3, m4] = get_means(max_muR, min_muR)
    m = min_muR + (max_muR-min_muR).*(rand(size(max_muR)));
    A = m;
    m1 = A(1,:);
    m2 = A(2,:);
    m3 = A(3,:);
    m4 = A(4,:);
end

function [sigma1, sigma2, sigma3, sigma4, theta] = get_covariance(alphaR,betaR)
    A = zeros(2,2,4);
    for i = 1:4 
        a = alphaR(1) + (alphaR(2)-alphaR(1))*rand;
        b = betaR(1) + (betaR(2)-betaR(1))*rand;
        theta(i) = 360.*rand;
        u = [cosd(theta(i)) sind(theta(i)); -sind(theta(i)) cosd(theta(i))];
        covar = u'*[a 0; 0 b]*u;
        A(:,:,i) = covar;
    end
    sigma1 = A(:,:,1);
    sigma2 = A(:,:,2);
    sigma3 = A(:,:,3);
    sigma4 = A(:,:,4);
end

function [f1,f2,f3,f4] = get_mvnpdf(mu1,mu2,mu3,mu4,Sigma1,Sigma2,Sigma3,Sigma4,x1,x2)
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

function [amp1, amp2, amp3, amp4] = get_amplitudes(delta, side)
    if side == 'L' % LeftHeavy
        amp1 = (0.5 + (0.2*rand))*delta(1); % amp1 lies within 50-70% of delta1
        amp2 = delta(1) - amp1;
        while 1
            amp3 = rand*delta(2);
            amp4 = delta(2) - amp3;
            % Check for the inequalities
            if(amp1 + amp4 > amp3 + amp2)
                break
            end
        end
    else
        amp3 = (0.5 + (0.2*rand))*delta(1); % amp1 lies within 50-70% of delta1
        amp4 = delta(1) - amp3;
        while 1
            amp1 = rand*delta(2);
            amp2 = delta(2) - amp1;
            % Check for the inequalities
            if(amp3 + amp2 > amp1 + amp4)
                break
            end
        end
    end
end

function Fp2 = randomize_halfM(Fp,r,n)
    % Comment below till 190 to not perform averaging or set flag = 1
    if r == 0
        for i = 1:n
            Fp = getNewMatrix(Fp);
        end
    end
    szFp = size(Fp,1);
   
    % Divide into half and shuffle
    % Comment below till 203 to not shuffle 
    % Fp2 = Fp
    p = fix(szFp/2);
    A = reshape(Fp(:,1:p),[],1);
    B = reshape(Fp(:,p+1:szFp),[],1);
   % %rng(1);
   randomIdxsA = randperm(numel(A));
   A_new = reshape(A(randomIdxsA),szFp,p);
   randomIdxsB = randperm(numel(B));
   B_new = reshape(B(randomIdxsB),szFp,szFp-p);
   Fp2 =  [A_new B_new];
end

function Fp2 = getNewMatrix(Fp)
    % Average of neighbouring columns and rows - (ODD issue)
     szFp = size(Fp,1); % Square matrix of odd size
    if mod(szFp,2) == 0
        iter = szFp/2;
    else
        iter = fix(szFp/2)+1;
    end
    Fp2 = zeros(iter);
    k = 1;
    for j = 1:szFp
        k=1;
        for i = 1:2:szFp
            if(i~=szFp)
                Fp2(j,k)= (Fp(j,i)+Fp(j,i+1))/2;
                k = k+1;
            else
                Fp2(j,k)=Fp(j,i);
            end
        end
    end
    for j = 1:iter % for all col
        for i = 1:2:szFp
            k = 1;
            if(i~=szFp)
                Fp2(k,j)= (Fp(i,j)+Fp(i+1,j))/2;
                k = k+1;
            else
                Fp2(k,j)=Fp(i,j);
            end
        end
    end
    Fp2(iter+1:szFp,:)=[];
 end