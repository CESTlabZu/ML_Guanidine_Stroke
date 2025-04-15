%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation code for tissue mimicking data. Please uncomment required 
% lines of text. Please run randidx line to generate a few random indices
% out of all possible simulated samples. Since the number of samples are
% very large, saving only few random samples can help save space. 
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
tic

% tt mean saturation power (uT)
tt_7pT =[0.9,0.95,1,1.05, 1.1]; % B1 shifts

% pulse parameters
pulseduration=5;
gauss=100;

% pool parameters
offppm= 300; 
sep1_7pT=3.6*offppm;
sep2_7pT=3*offppm;
sep3_7pT=2*offppm
sep4_7pT=-1.6*offppm;
sep5_7pT=-3.3*offppm;


% relaxations
R1S=1/1.5;
R2S2=1/0.015;
R2S4=1/0.001;
R2S5=1/0.0005;
R1M=1/1.5;

% frequency offset
maxi=1500;
step=50;
offset= -maxi:step:maxi;
k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];

k_7pT = [k-40;k-20;k;k+20;k+40]'; % B0 shifts

satangle=tt_7pT*42.6*360*pulseduration;

% constant sample parameter

fs4=0.003;
ksw4=50;
ksw5=20;
kmw=25;


num_T1W=4;
num_T2W=3;
num_T2S1=3;
num_T2S3=4;
num_T2M=3;

num_fs1=3;
num_fs2=3;
num_fs3=3;
num_fs5=3;
num_fm=3;

num_ksw1=3;
num_ksw2=3;
num_ksw3=3;

num_tt = 5;
num_B0 = 5;

T1W_matrix=[1.6,1.8,2.0,2.2];
T2W_matrix=[60, 80, 100]*0.001;
T2S1_matrix = [0.002, 0.003, 0.004];
T2S3_matrix=[0.007, 0.009, 0.011, 0.013];
T2M_matrix = [0.00003, 0.00005,0.00007];

fs1_matrix = [0.0005,0.0009,0.0013];
fs2_matrix = [0.6,0.8,1.0]*0.003;
fs3_matrix =[0.0003, 0.0004, 0.0005];
fs5_matrix = [0.002,0.008,0.014];
fm_matrix=[0.05,0.08, 0.11];

ksw1_matrix = [50, 80, 110];
ksw2_matrix = [3000,5000,7000];
ksw3_matrix = [300, 500, 700];

i=1;


for ii_T1W=1:num_T1W
    ii_T1W
    R1W=1./T1W_matrix(ii_T1W);
    for ii_T2W=1:num_T2W
        ii_T2W
        R2W=1./T2W_matrix(ii_T2W);
        for ii_T2S1=1:num_T2S1
            R2S1 = 1./T2S1_matrix(ii_T2S1);
            for ii_T2S3=1:num_T2S3
                R2S3 = 1./T2S3_matrix(ii_T2S3);
                for ii_T2M = 1:num_T2M
                    R2M = 1./T2M_matrix(ii_T2M);
                    for ii_fs1=1:num_fs1
                        fs1=fs1_matrix(ii_fs1);
                        for ii_fs2=1:num_fs2
                            fs2=fs2_matrix(ii_fs2);
                            for ii_fs3=1:num_fs3
                                fs3=fs3_matrix(ii_fs3);
                                for ii_fs5=1:num_fs5
                                    fs5=fs5_matrix(ii_fs5);
                                    for ii_fm=1:num_fm
                                        fm=fm_matrix(ii_fm);
                                        for ii_ksw1=1:num_ksw1
                                            ksw1=ksw1_matrix(ii_ksw1);
                                            for ii_ksw2=1:num_ksw2
                                                ksw2=ksw2_matrix(ii_ksw2);
                                                for ii_ksw3=1:num_ksw3
                                                    ksw3=ksw3_matrix(ii_ksw3);
                                                    for ii_tt=1:num_tt
                                                        tt_shift = satangle(ii_tt); 
                                                        for ii_kk = 1:num_B0
                                                            B0_shift = k_7pT(:,ii_kk); 
                                                  
                                 
R1W_cal_obs=(R1W+(fm*R1M))./(1+fm); 
R1W_cal_matrix(ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk)=R1W_cal_obs;                                

a25mspulse = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, tt_shift, 1, 2, 1, .00, 1, 1, B0_shift*2*pi, 1);
a25mspulse_ref = runsteadysimgauss(ksw1,ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, 0,fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, tt_shift, 1, 2, 1, .00, 1, 1, B0_shift*2*pi, 1);

% without B0 and B1 shift
a25mspulse_ns = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, satangle(3), 1, 2, 1, .00, 1, 1, k_7pT(:,3)*2*pi, 1);
a25mspulse_ref_ns = runsteadysimgauss(ksw1,ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, 0,fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, satangle(3), 1, 2, 1, .00, 1, 1, k_7pT(:,3)*2*pi, 1);

aa_7pT(:,ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk)=a25mspulse(:,6);
aa_7pT_ref(:,ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk)=a25mspulse_ref(:,6);

Slab = a25mspulse(:,6);
Sref = a25mspulse_ref(:,6);

Slab_ns = a25mspulse_ns(:,6);
Sref_ns = a25mspulse_ref_ns(:,6);

S0 = 1;

fm_cal_matrix(ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk)=fm;
params_matrix(:, ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk) = [R1W_cal_obs; R2W; R2S1; R2S3; R2M; fs1; fs2; fs3; fs5; fm; ksw1; ksw2; ksw3; tt_7pT(ii_tt); k_7pT(35,ii_kk)];
 
mtr = reshape(((1-Slab_ns)-(1-Sref_ns)), [69 1]);
mtr_amp(ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk) = mtr(23,:);
mtr_width(ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk)= fwhm2(mtr,k_7pT(:,3));

arex=((1./Slab_ns) - (1./Sref_ns)).*R1W_cal_obs*(1+fm);
arex_amp(ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk) = arex(23);
arex_width(ii_T1W,ii_T2W,ii_T2S1, ii_T2S3, ii_T2M, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm, ii_ksw1, ii_ksw2, ii_ksw3, ii_tt,ii_kk)= fwhm2(arex(1:30),k_7pT(:,3));
 
X = sprintf("-------------------------------------%d",i);
disp(X)
i = i+1;
                                          end
                                       end
                                    end
                                 end
                              end
                           end
                        end
                     end
                  end
               end
            end
         end
      end
   end
end



matrix_output1(:,1)=reshape(mtr_amp,  [num_T1W*num_T2W*num_T2S1*num_T2S3*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3*num_tt*num_B0 1]);     
matrix_output1(:,2)=reshape(mtr_width,  [num_T1W*num_T2W*num_T2S1*num_T2S3*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3*num_tt*num_B0  1]);

matrix_output2(:,1)=reshape(arex_amp,  [num_T1W*num_T2W*num_T2S1*num_T2S3*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3*num_tt*num_B0  1]);     
matrix_output2(:,2)=reshape(arex_width,  [num_T1W*num_T2W*num_T2S1*num_T2S3*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3*num_tt*num_B0  1]);

matrix_input_all=reshape(aa_7pT(:,:),  [69, num_T1W*num_T2W*num_T2S1*num_T2S3*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3*num_tt*num_B0]);   
R1W_cal_matrix_output_all=reshape(R1W_cal_matrix,  [num_T1W*num_T2W*num_T2S1*num_T2S3*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3*num_tt*num_B0  1]);
fm_cal_matrix_output_all=reshape(fm_cal_matrix,  [num_T1W*num_T2W*num_T2S1*num_T2S3*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3*num_tt*num_B0  1]);

tm_randidx = randperm(size(matrix_input_all,2), 5000); % saving only few indices can help saving space

toc