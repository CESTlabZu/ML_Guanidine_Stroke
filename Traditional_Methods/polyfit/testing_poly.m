%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code for polynomial fitting. Please uncomment required lines of text.
% 
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please load your input CEST files, observed R1 and fm map.
load('sample_input.mat');
load('sample_r1.mat');

% Creating frequency offset
max=1500;
step=50;
offset= -max:step:max;
k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]/300;

FitParam.PeakOffset = -2;
outputaptspec1 = zeros([1 1 1 69]);

for i=1:64
    for j=1:64
        for kk = 1:4
            if(squeeze(sample_input(i,j,1,kk)) == 0)
                outputaptspec1(i,j,kk,:) = zeros([1,1,1,69]);
            else
                FitParam.R1 = r1(i,j,kk);                               
                FitParam.fm_m = fm(i,j,kk); 
                Zspectra_INPUT = flip(squeeze(sample_input(i,j,:,kk)));
                fitresult = poly(k,Zspectra_INPUT', FitParam); % change to poly_opt if optimized frequency 
                                                               % offset fitting wanted
                background = fitresult.Background;
                sref = background';
                outputaptspec1(i,j,kk,:) = fitresult.Arex;  
                % outputaptspec1(i,j,kk,17:26) = fitresult.Arex;   % please uncomment this line for optimized frequency 
                                                                   % offset fitting 
            end
        end
    end
    i
end

% Plotting fitted results
imagesc(outputaptspec1(:,:,3,23))
clim([0 0.06]) % please adjust limits for optimized offsets after polynomial fitting~!