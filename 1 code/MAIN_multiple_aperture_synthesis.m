%--------------------------------------------Project Information----------------------------------------------%
% Project Name: Depth-of-focus extension in optical coherence tomography via multiple aperture synthesis
% Institution: Nanyang Technological University
% Creator: Bo En
% Creation Date: 31-March-2017
% Description: This algorithm aims to digitally synthesize five distinctive cross-sectional images of 6um beads
%              captured via five distinctive apertures into one brand new DOF-extended image.
% Reference paper: E. Bo, Y. Luo, S. Chen, X. Liu, N. Wang, X. Ge, X. Wang, S. Chen, S. Chen, and J. Li,
%                  "Depth-of-focus extension in optical coherence tomography via multiple aperture synthesis,"
%                  Optica 4, 701-706 (2017). https://doi.org/10.1364/OPTICA.4.000701
% Modification history: code simplification for publication online on 27-Jan-2018.
%--------------------------------------------------Main Code--------------------------------------------------%
clear all; close all; clc; warning off;
addpath('C:\Users\boen0001\Desktop\Depth-of-focus extension in OCT via multiple aperture synthesis\1 code');
%% parameter definition
c = 1024;
pn = 4096;
NFFT = 4096*1;
z_st = 76; z_end = 695; x_st = 274; x_end = 785;
z_sft_st = -2; z_sft_end = 2; pha_num = 4;

%% capture 20 interference fringes, run dispersion_compensation.m
dispersion_compensation;
dir_set(2);
load('imfirst');
xip = CalStru.xi;
MA = CalStru.MAmean;
CArray = repmat(CalStru.CArray,1,c);

%% Read data
dir_set(3);
load('A1'); load('A2'); load('A3'); load('A4'); load('A5');

%% Bgn process
[ref1] = bgn_process(A1,4095,c);
[ref2] = bgn_process(A2,4095,c);
[ref3] = bgn_process(A3,4095,c);
[ref4] = bgn_process(A4,4095,c);
[ref5] = bgn_process(A5,4095,c);

%% Output original images
B1 = interp1(MA,(A1-ref1),xip);
B2 = interp1(MA,(A2-ref2),xip);
B3 = interp1(MA,(A3-ref3),xip);
B4 = interp1(MA,(A4-ref4),xip);
B5 = interp1(MA,(A5-ref5),xip);
BSCAN1 = single(abs(fft(B1.*CArray,NFFT))); BSCAN1 = BSCAN1(z_st:z_end,x_st:x_end);
BSCAN2 = single(abs(fft(B2.*CArray,NFFT))); BSCAN2 = BSCAN2(z_st:z_end,x_st:x_end);
BSCAN3 = single(abs(fft(B3.*CArray,NFFT))); BSCAN3 = BSCAN3(z_st:z_end,x_st:x_end);
BSCAN4 = single(abs(fft(B4.*CArray,NFFT))); BSCAN4 = BSCAN4(z_st:z_end,x_st:x_end);
BSCAN5 = single(abs(fft(B5.*CArray,NFFT))); BSCAN5 = BSCAN5(z_st:z_end,x_st:x_end);
dir_set(3);
print_image('BSCAN1.bmp',BSCAN1);
print_image('BSCAN2.bmp',BSCAN2);
print_image('BSCAN3.bmp',BSCAN3);
print_image('BSCAN4.bmp',BSCAN4);
print_image('BSCAN5.bmp',BSCAN5);

%% MAS STEP 1: axial shift
[z_sft1, z_sft2, z_sft3, z_sft4, z_sft5] = opt_z_sft(B1, B2, B3, B4, B5, z_sft_st, z_sft_end, z_st, z_end, x_st, x_end);
B1_z_sfted = fsft(B1,z_sft1,0,CArray,NFFT);
B2_z_sfted = fsft(B2,z_sft2,0,CArray,NFFT);
B3_z_sfted = fsft(B3,z_sft3,0,CArray,NFFT);
B4_z_sfted = fsft(B4,z_sft4,0,CArray,NFFT);
B5_z_sfted = fsft(B5,z_sft5,0,CArray,NFFT);
STEP1 = abs(B1_z_sfted(z_st:z_end,x_st:x_end) + B2_z_sfted(z_st:z_end,x_st:x_end)...
    + B3_z_sfted(z_st:z_end,x_st:x_end) + B4_z_sfted(z_st:z_end,x_st:x_end)...
    + B5_z_sfted(z_st:z_end,x_st:x_end))/5;
dir_set(4);
print_image('MAS_STEP1.bmp',STEP1);

%% MAS STEP 2: correct phase
pha1 = zeros(z_end - z_st + 1,1); pha2 = pha1; pha3 = pha1; pha4 = pha1; pha5 = pha1;
for depth = z_st:z_end
    depth;
    [opt_pha1, opt_pha2, opt_pha3, opt_pha4, opt_pha5] = opt_pha(B1_z_sfted(depth,x_st:x_end), B2_z_sfted(depth,x_st:x_end), B3_z_sfted(depth,x_st:x_end), B4_z_sfted(depth,x_st:x_end), B5_z_sfted(depth,x_st:x_end), pha_num);
    pha1(depth-z_st+1,1) = opt_pha1; pha2(depth-z_st+1,1) = opt_pha2; pha3(depth-z_st+1,1) = opt_pha3; pha4(depth-z_st+1,1) = opt_pha4; pha5(depth-z_st+1,1) = opt_pha5;
end
pha1 = repmat(pha1,1,x_end-x_st+1); pha2 = repmat(pha2,1,x_end-x_st+1);pha3 = repmat(pha3,1,x_end-x_st+1);pha4 = repmat(pha4,1,x_end-x_st+1);pha5 = repmat(pha5,1,x_end-x_st+1);
B1_coted = B1_z_sfted(z_st:z_end,x_st:x_end) .* exp(-1i * pha1);
B2_coted = B2_z_sfted(z_st:z_end,x_st:x_end) .* exp(-1i * pha2);
B3_coted = B3_z_sfted(z_st:z_end,x_st:x_end) .* exp(-1i * pha3);
B4_coted = B4_z_sfted(z_st:z_end,x_st:x_end) .* exp(-1i * pha4);
B5_coted = B5_z_sfted(z_st:z_end,x_st:x_end) .* exp(-1i * pha5);
STEP2 = abs(B1_coted + B2_coted + B3_coted + B4_coted + B5_coted)/5;
print_image('MAS_STEP2.bmp',STEP2);