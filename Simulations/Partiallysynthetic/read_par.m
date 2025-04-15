clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation code for partially synthetic data. Please uncomment required 
% lines of text. 
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load 6-pool fitted amines/guan and MT lineshapes
% load('...');
% arex_amine = ... % shape [69 1]
% arex_MT = ... % shape [69 1]

% tt mean saturation power (uT)
tt_7pT =[0.9,0.95,1,1.05, 1.1]; % B1 shifts


% pool paramters
offppm = 300;
sep1=3.6*offppm;
sep2=3*offppm;
sep3=2*offppm
sep4=-1.6*offppm;
sep5=-3.3*offppm;

% mfit-poly to obtain amine lineshape - measured amine component
maxi=1500;
step=50;
offset= -maxi:step:maxi;
k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
FitParam.PeakOffset = -3;

% smooth curve to denoise lineshape
Zdenoised_normal = mean(arex_amine_n,2);
Zdenoised_normal = abs(Zdenoised_normal);
Zdenoised_normal(1:5) = smoothdata(Zdenoised_normal(1:5),'gaussian', 5);
Zdenoised_normal(6:31) = smoothdata(Zdenoised_normal(6:31),'sgolay',5);
Zdenoised_normal(41:69) = smoothdata(Zdenoised_normal(41:69),'gaussian',25);

fitresult = mfit_poly(k,Zdenoised_normal', FitParam);
background_normal = fitresult.Background;
amine_normal = [Zdenoised_normal(1:16)', background_normal, Zdenoised_normal(32:69)'];


Zdenoised_lesion = mean(arex_amine_l,2);
Zdenoised_lesion = abs(Zdenoised_lesion);
Zdenoised_lesion(1:8) = smoothdata(Zdenoised_lesion(1:8),'gaussian', 5);
Zdenoised_lesion(9:31) = smoothdata(Zdenoised_lesion(9:31),'sgolay',5);
Zdenoised_lesion(40:69) = smoothdata(Zdenoised_lesion(40:69),'gaussian',25);

fitresult = mfit_poly(k,Zdenoised_lesion', FitParam);
background_lesion = fitresult.Background;
amine_lesion = [Zdenoised_lesion(1:16)', background_lesion, Zdenoised_lesion(32:69)'];

% measured MT component
m_MT_n = mean(arex_MT_n,2);
m_MT_l = mean(arex_MT_l,2);

% constant sample parameters for simulations
R2S1=1/0.002;
R2S5=1/0.0005;
ksw5 = 20;
R1M = 1/1.5;

% simulated components
fm_n=mean(fm_n,2);
R1W_n=mean(r1w_n,2);
fm_l=mean(fm_l,2);
R1W_l=mean(r1w_l,2);

num_T1W=3;
num_T2W=3;
num_T2S3=3;

num_fs1=3;
num_fs2=5;
num_fs3=4;
num_fs5=3;
num_fm=3;

num_ksw1=5;
num_ksw3=3;

num_tt = 5;
num_B0 = 5;

k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];

k_7pT = [k-40;k-20;k;k+20;k+40];

T1W_matrix=[1.5,2.0,2.5];
T2W_matrix=[50, 80, 110]*0.001;
T2S3_matrix=[0.005, 0.010, 0.015];
 
fs1_matrix=[0.0006,0.0010,0.0014];
fs2_matrix=[0.5,0.75,1,1.25,1.5];
fs3_matrix=[0.0001,0.0003,0.0005,0.0007];
fs5_matrix=[0.002,0.008,0.014];
fm_matrix=[0.8, 1.0, 1.2];

ksw1_matrix = [40,70,100,130,160];
ksw3_matrix = [300, 500, 700];



i=1;
 
for ii_T1W=1:num_T1W
    ii_T1W
    R1W_cal=1./T1W_matrix(ii_T1W);
    for ii_T2W=1:num_T2W
        ii_T2W
        R2W_cal=1./T2W_matrix(ii_T2W);
        for ii_T2S=1:num_T2S3
            R2S3_cal = 1./T2S3_matrix(ii_T2S);
            for ii_fs1=1:num_fs1
                fs1_cal=fs1_matrix(ii_fs1);
                for ii_fs2=1:num_fs2
                    fs2_cal=fs2_matrix(ii_fs2);
                    for ii_fs3 = 1:num_fs3
                        fs3_cal =fs3_matrix(ii_fs3);
                        for ii_fs5=1:num_fs5
                            fs5_cal=fs5_matrix(ii_fs5);
                            for ii_fm=1:num_fm
                                fm_cal=fm_matrix(ii_fm);
                                for ii_ksw1 = 1:num_ksw1
                                    ksw1 = ksw1_matrix(ii_ksw1);
                                    for ii_ksw3 = 1:num_ksw3
                                        ksw3_cal = ksw3_matrix(ii_ksw3);
                                        for ii_tt = 1:num_tt
                                            tt_shift = tt_7pT(ii_tt); 
                                            for ii_kk = 1:num_B0
                                                B0_shift = k_7pT(ii_kk,:);
                               


% simulations with normal components
R1W_cal_obs_n=(R1W_cal+fm_cal*fm_n*R1M)./(1+fm_cal*fm_n);
cal_Lorentzian1_cal=(fs1_cal.*ksw1.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S1+ksw1)*ksw1+ksw1./(R2S1+ksw1).*((B0_shift+sep1)*2*pi).^2));
cal_Lorentzian2_n_cal=fs2_cal.*(interp1(k_7pT(3,:),amine_normal,B0_shift))*((tt_shift)^2);
cal_Lorentzian3_n_cal = (fs3_cal.*ksw3_cal.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S3_cal+ksw3_cal)*ksw3_cal+ksw3_cal./(R2S3_cal+ksw3_cal).*((B0_shift+sep3)*2*pi).^2));
cal_Lorentzian5_n_cal = (fs5_cal.*ksw5.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((B0_shift+sep5)*2*pi).^2));
cal_Lorentzian6_n_cal=fm_cal*(interp1(k_7pT(3,:),m_MT_n,B0_shift))*((tt_shift)^2);
cal_eff_cal=R1W_cal_obs_n.*((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2)+R2W_cal.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2);
sscal_n = R1W_cal_obs_n./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm_n)+cal_Lorentzian2_n_cal./(1+fm_cal*fm_n)+cal_Lorentzian3_n_cal./(1+fm_cal*fm_n)+cal_Lorentzian5_n_cal./(1+fm_cal*fm_n)+cal_Lorentzian6_n_cal).*(((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2));

SS_cal_n(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm, ii_ksw1, ii_ksw3, ii_tt, ii_kk)=sscal_n;
SS_cal_value_ref_n=R1W_cal_obs_n./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm_n)+cal_Lorentzian2_n_cal./(1+fm_cal*fm_n)+0./(1+fm_cal*fm_n)+cal_Lorentzian5_n_cal./(1+fm_cal*fm_n)+cal_Lorentzian6_n_cal).*(((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2));

% without B0 and B1 shifts
cal_Lorentzian1_cal_ns=(fs1_cal.*ksw1.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+(R2S1+ksw1)*ksw1+ksw1./(R2S1+ksw1).*((k+sep1)*2*pi).^2));
cal_Lorentzian2_n_cal_ns=amine_normal;
cal_Lorentzian3_n_cal_ns = (fs3_cal.*ksw3_cal.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+(R2S3_cal+ksw3_cal)*ksw3_cal+ksw3_cal./(R2S3_cal+ksw3_cal).*((k+sep3)*2*pi).^2));
cal_Lorentzian5_n_cal_ns = (fs5_cal.*ksw5.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((k+sep5)*2*pi).^2));
cal_Lorentzian6_n_cal_ns= m_MT_n';
cal_eff_cal_ns=R1W_cal_obs_n.*((k)*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2)+R2W_cal.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2);
sscal_n_ns = R1W_cal_obs_n./(cal_eff_cal_ns+cal_Lorentzian1_cal_ns./(1+fm_cal*fm_n)+cal_Lorentzian2_n_cal_ns./(1+fm_cal*fm_n)+cal_Lorentzian3_n_cal_ns./(1+fm_cal*fm_n)+cal_Lorentzian5_n_cal_ns./(1+fm_cal*fm_n)+cal_Lorentzian6_n_cal_ns).*(((k)*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2));
SS_cal_n_ns(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm, ii_ksw1, ii_ksw3, ii_tt, ii_kk)=sscal_n_ns;
SS_cal_value_ref_n_ns=R1W_cal_obs_n./(cal_eff_cal_ns+cal_Lorentzian1_cal_ns./(1+fm_cal*fm_n)+cal_Lorentzian2_n_cal_ns./(1+fm_cal*fm_n)+0./(1+fm_cal*fm_n)+cal_Lorentzian5_n_cal_ns./(1+fm_cal*fm_n)+cal_Lorentzian6_n_cal_ns).*(((k)*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2));

% simulation with lesion components
R1W_cal_obs_l=(R1W_cal+fm_cal*fm_l*R1M)./(1+fm_cal*fm_l);
cal_Lorentzian1_cal=(fs1_cal.*ksw1.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S1+ksw1)*ksw1+ksw1./(R2S1+ksw1).*((B0_shift+sep1)*2*pi).^2));
cal_Lorentzian2_l_cal=fs2_cal.*(interp1(k_7pT(3,:),amine_lesion,B0_shift))*((tt_shift)^2);
cal_Lorentzian3_l_cal = (fs3_cal.*ksw3_cal.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S3_cal+ksw3_cal)*ksw3_cal+ksw3_cal./(R2S3_cal+ksw3_cal).*((B0_shift+sep3)*2*pi).^2));
cal_Lorentzian5_l_cal = (fs5_cal.*ksw5.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((B0_shift+sep5)*2*pi).^2));
cal_Lorentzian6_l_cal=fm_cal*(interp1(k_7pT(3,:),m_MT_l,B0_shift))*((tt_shift)^2);
cal_eff_cal=R1W_cal_obs_l.*((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2)+R2W_cal.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2);
sscal_l = R1W_cal_obs_l./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm_l)+cal_Lorentzian2_l_cal./(1+fm_cal*fm_l)+cal_Lorentzian3_l_cal./(1+fm_cal*fm_l)+cal_Lorentzian5_l_cal./(1+fm_cal*fm_l)+cal_Lorentzian6_l_cal).*(((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2));
SS_cal_l(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm, ii_ksw1, ii_ksw3, ii_tt, ii_kk)=sscal_l;
SS_cal_value_ref_l=R1W_cal_obs_l./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm_l)+cal_Lorentzian2_l_cal./(1+fm_cal*fm_l)+0./(1+fm_cal*fm_l)+cal_Lorentzian5_l_cal./(1+fm_cal*fm_l)+cal_Lorentzian6_l_cal).*(((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2));

% without B0 and B1 shifts
cal_Lorentzian1_cal_ns=(fs1_cal.*ksw1.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+(R2S1+ksw1)*ksw1+ksw1./(R2S1+ksw1).*((k+sep1)*2*pi).^2));
cal_Lorentzian2_l_cal_ns=amine_lesion;
cal_Lorentzian3_l_cal_ns = (fs3_cal.*ksw3_cal.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+(R2S3_cal+ksw3_cal)*ksw3_cal+ksw3_cal./(R2S3_cal+ksw3_cal).*((k+sep3)*2*pi).^2));
cal_Lorentzian5_l_cal_ns = (fs5_cal.*ksw5.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((k+sep5)*2*pi).^2));
cal_Lorentzian6_l_cal_ns= m_MT_l';
cal_eff_cal_ns=R1W_cal_obs_l.*((k)*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2)+R2W_cal.*(tt_7pT(3).*42.6*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2);
sscal_l_ns = R1W_cal_obs_l./(cal_eff_cal_ns+cal_Lorentzian1_cal_ns./(1+fm_cal*fm_l)+cal_Lorentzian2_l_cal_ns./(1+fm_cal*fm_l)+cal_Lorentzian3_l_cal_ns./(1+fm_cal*fm_l)+cal_Lorentzian5_l_cal_ns./(1+fm_cal*fm_l)+cal_Lorentzian6_l_cal_ns).*(((k)*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2));
SS_cal_l_ns(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm, ii_ksw1, ii_ksw3, ii_tt, ii_kk)=sscal_l_ns;
SS_cal_value_ref_l_ns=R1W_cal_obs_l./(cal_eff_cal_ns+cal_Lorentzian1_cal_ns./(1+fm_cal*fm_l)+cal_Lorentzian2_l_cal_ns./(1+fm_cal*fm_l)+0./(1+fm_cal*fm_l)+cal_Lorentzian5_l_cal_ns./(1+fm_cal*fm_l)+cal_Lorentzian6_l_cal_ns).*(((k)*2*pi).^2./((tt_7pT(3).*42.6*2*pi).^2+((k)*2*pi).^2));



R1W_cal_matrix_n(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)=R1W_cal_obs_n;
fm_cal_matrix_n(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)=fm_cal*fm_n;

R1W_cal_matrix_l(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)=R1W_cal_obs_l;
fm_cal_matrix_l(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)=fm_cal*fm_l;

params_matrix(:,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk) = [R1W_cal_obs_n; R2W_cal;R2S3_cal; fs1_cal; fs2_cal; fs3_cal; fs5_cal; fm_cal; ksw3_cal; tt_shift; B0_shift(35)];

mtr_n = SS_cal_value_ref_n_ns -  sscal_n_ns;
mtr_amp_n(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk) = mtr_n(1,23);
mtr_width_n(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)= fwhm2(mtr_n,k);

arex_n=((1./sscal_n_ns) - (1./SS_cal_value_ref_n_ns)).*(R1W_cal_obs_n).*(1+(fm_cal*fm_n));
arex_amp_n(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk) = arex_n(1,23);
arex_width_n(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)= fwhm2(arex_n,k);

mtr_l = SS_cal_value_ref_l_ns -  sscal_l_ns;
mtr_amp_l(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk) = mtr_l(1,23);
mtr_width_l(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)= fwhm2(mtr_l,k);

arex_l=((1./sscal_l_ns) - (1./SS_cal_value_ref_l_ns)).*(R1W_cal_obs_l).*(1+(fm_cal*fm_l));
arex_amp_l(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk) = arex_l(1,23);
arex_width_l(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5, ii_fm,ii_ksw3, ii_tt, ii_kk)= fwhm2(arex_l,k);

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

matrix_input_all_n(:,:)=reshape(SS_cal_n,  [length(k') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0]); 
R1W_cal_matrix_output_all_n=reshape(R1W_cal_matrix_n,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);
fm_cal_matrix_output_all_n=reshape(fm_cal_matrix_n,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);
params_matrix = reshape(params_matrix, [11 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0]);
matrix_MTR_output_n(:,1)=reshape(mtr_amp_n,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);   
matrix_MTR_output_n(:,2)=reshape(mtr_width_n,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);
matrix_AREX_output_n(:,1)=reshape(arex_amp_n,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);   
matrix_AREX_output_n(:,2)=reshape(arex_width_n,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);


matrix_input_all_l(:,:)=reshape(SS_cal_l,  [length(k') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0]); 
R1W_cal_matrix_output_all_l=reshape(R1W_cal_matrix_l,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);
fm_cal_matrix_output_all_l=reshape(fm_cal_matrix_l,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);
matrix_MTR_output_l(:,1)=reshape(mtr_amp_l,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);   
matrix_MTR_output_l(:,2)=reshape(mtr_width_l,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);
matrix_AREX_output_l(:,1)=reshape(arex_amp_l,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);   
matrix_AREX_output_l(:,2)=reshape(arex_width_l,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw3*num_tt*num_B0 1]);

