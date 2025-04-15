%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code for 7-pool Lorentzian fitting. Please uncomment required lines 
% of text.
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% load CEST Z spectrum, R1W and fm values from simulations or measured
% data.
load('sample_input.mat');
load('sample_r1.mat');

index = [1, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,...
        28, 31, 37, 40, 53, 63]+1; % optimized offset indices (python + 1)

% required initial parameters
max=1500;
step=50;
offset= -max:step:max;
k_7pT=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
k_7pT=k_7pT';

for ii=1:64
    for jj=1:64
        for kk=1:4
            if(sample_input(ii,jj,1,kk) == 0)
                mor_AREX_amine(ii,jj,kk,:) = zeros(1,69);
                mor_AREX_MT(ii,jj,kk,:) = zeros(1,69);
                mor_AREX_guan(ii,jj,kk,:)  = zeros(1,69);
                mor_AREX_amide(ii,jj,kk,:) = zeros(1,69);
                mor_AREX_NOE3p5(ii,jj,kk,:) = zeros(1,69);
                mor_AREX_NOE1p6(ii,jj,kk,:) = zeros(1,69);
            else
                matrix_input = flip(squeeze(sample_input(ii,jj,:,kk)));

                % uncomment for fitting with optimized offsets
                % matrix_input = matrix_input(index);

                R1W = squeeze(r1(ii,jj,kk));
                fm_m = squeeze(fm(ii,jj,kk));
    
                sig=(1-matrix_input); 
                R1W_AREX=R1W; % loaded R1W values
                fm_AREX=fm_m; % loaded PSR values
               
                x=k_7pT;
                % x = k_7pT(index); % please uncomment for fitting with
                                    % optimized frequency offset
                % beta = [A, Off, W]
                % beta = [Water, amide, amine, guanidine, NOE-1.6, NOE-3.5, MT]
                beta0= [0.9, 0, 420,           0.025, -1050, 150,       0.1, -900, 450,   0.01,-600,300,       0.001, 450, 300,         0.02, 1050, 900,       0.1, 0, 7500]; % initial test
                lb=[  0.02, -300, 30,          0, -1200,120,       0, -750, 150,          0, -750, 150,      0, 300, 0,              0, 750, 300,           0, -1200, 3000]; % lower bound
                ub=[ 1, 300,   3000,      0.2, -900, 900,         0.3,-1050, 1500,        0.1,-450,600,      0.2, 600, 450,       1, 1350, 1500,         1, 1200, 30000]; % upper bound
                
                Delta=[1]; 
                options=optimset('lsqcurvefit') ; 
                options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;
                
                [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                    lsqcurvefit(@matsolv_7pool, beta0, x, sig, lb, ub, options, Delta) ;
                
                
                x = k_7pT;
                % amide
                beta_amide=beta;
                sig_simur_amide=matsolv_7pool(beta_amide,x,Delta);
                beta_amide(4)=0;
                sig_simur_ref_amide=matsolv_7pool(beta_amide,x,Delta);

                % save required values!
                mor_AREX_amide(ii,jj,kk,:)=(1./(1-sig_simur_amide)-1./(1-sig_simur_ref_amide))*R1W_AREX*(1+fm_AREX);

                % amine
                beta_amine=beta;
                sig_simur_amine=matsolv_7pool(beta_amine,k_7pT,Delta);
                beta_amine(7)=0;
                sig_simur_ref_amine=matsolv_7pool(beta_amine,k_7pT,Delta);

                mor_AREX_amine(ii,jj,kk,:)=(1./(1-sig_simur_amine)-1./(1-sig_simur_ref_amine))*R1W_AREX*(1+fm_AREX);


                % guanidine
                beta_guan=beta;
                sig_simur_guan=matsolv_7pool(beta_guan,k_7pT,Delta);
                beta_guan(10)=0;

                sig_simur_ref_guan=matsolv_7pool(beta_guan,k_7pT,Delta);
                mor_AREX_guan(ii,jj,kk,:)=(1./(1-sig_simur_guan)-1./(1-sig_simur_ref_guan))*R1W_AREX*(1+fm_AREX);

                % MT
                beta_MT=beta;
                sig_simur_MT=matsolv_7pool(beta_MT,x,Delta);
                beta_MT(19)=0;
                sig_simur_ref_MT=matsolv_7pool(beta_MT,x,Delta);

                mor_MT=(sig_simur_MT-sig_simur_ref_MT);
                mor_AREX_MT(ii,jj,kk,:) = (mor_MT./(1-mor_MT))*R1W_AREX;

                % NOE3p5
                beta_NOE3p5=beta;
                sig_simur_NOE3p5=matsolv_7pool(beta_NOE3p5,x,Delta);
                beta_NOE3p5(16)=0;
                sig_simur_ref_NOE3p5=matsolv_7pool(beta_NOE3p5,x,Delta);

                mor_AREX_NOE3p5(ii,jj,kk,:)=(1./(1-sig_simur_NOE3p5)-1./(1-sig_simur_ref_NOE3p5))*R1W_AREX*(1+fm_AREX);

                % NOE1p6
                beta_NOE1p6=beta;
                sig_simur_NOE1p6=matsolv_7pool(beta_NOE1p6,x,Delta);
                beta_NOE1p6(13)=0;
                sig_simur_ref_NOE1p6=matsolv_7pool(beta_NOE1p6,x,Delta);

                mor_AREX_NOE1p6(ii,jj,kk,:)=(1./(1-sig_simur_NOE1p6)-1./(1-sig_simur_ref_NOE1p6))*R1W_AREX*(1+fm_AREX);
            end
            
        end
    end
    sprintf("----------------------- %d",ii)
end

% plotting sample fitting!
imagesc(mor_AREX_guan(:,:,2,23))
clim([0 0.06])