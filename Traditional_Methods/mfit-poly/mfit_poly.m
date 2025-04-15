%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code for mfit-poly fitting. Please uncomment required lines of text.
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% load 6-pool Lorentzian fitted amine spectrum! For fitting using
% optimized offsets, please load 6-pool Lorentzian fitting using
% optimized offsets
load('sample_input.mat')


% required parameters
max=1500;
step=50;
offset= -max:step:max;

k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]/300;
FitParam.PeakOffset = -2;
outputaptspec1 = zeros([1 1 1 69]);

for i=1:64
    for j=1:64
        for kk = 1:4
            if(squeeze(sample_input(i,j,kk,1)) == 0)
                outputaptspec1(i,j,kk,:) = zeros([1,1,1,69]);
            else
                Zspectra_INPUT = squeeze(sample_input(i,j,kk,:));
                fitresult = mfitpoly(k,Zspectra_INPUT', FitParam);
                background = fitresult.Background';
                outputaptspec1(i,j,kk,15:26) = (Zspectra_INPUT(15:26) - background);
            end
        end
    end
    i
end

% plotting sample fitting
imagesc(outputaptspec1(:,:,2,23))
clim([-0.03 0.06])
