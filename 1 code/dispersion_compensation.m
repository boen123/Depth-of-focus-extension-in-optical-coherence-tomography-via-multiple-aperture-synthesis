function [] = dispersion_compensation()
% this function aims to perform dispersion compensation and spectrum estimation
% reference: X. Liu, S. Chen, D. Cui, X. Yu, and L. Liu, "Spectral estimation optical coherence tomography for axial super-resolution," Optics express 23, 26521-26532 (2015).

%% parameter defination
FFT_NUM = 32;
NFFT = 4096*FFT_NUM;
pn = 4096;             
c = 1024;              
NUM = 20;
BGN_NUM = NUM;
DISPERSION_NUM = NUM;
PLOT_NUM = NUM;
x_scale_rate = 0.5/223/FFT_NUM;

addpath('C:\Users\boen0001\Desktop\Depth-of-focus extension in OCT via multiple aperture synthesis\1 code');
cd('C:\Users\boen0001\Desktop\Depth-of-focus extension in OCT via multiple aperture synthesis\2 dispersion compensation');
%% Background process
sum = 0;
for index=1:BGN_NUM
    f = importdata(['f',num2str(index),'.mat']);
    sum = sum+f;
end
s = mean(sum/BGN_NUM,2);
r = 0;

%% Dispersion compensation
Ang_Cutting_Limit = [100,3450];
T1 = round((Ang_Cutting_Limit(2)-Ang_Cutting_Limit(1))/3);
T2 = 2*round((Ang_Cutting_Limit(2)-Ang_Cutting_Limit(1))/3);
parfor frame = 1:DISPERSION_NUM/2
    display(['Calculate MA :',num2str(frame)]);
    f1 = importdata(['f',num2str(frame),'.mat']);
    f2 = importdata(['f',num2str(DISPERSION_NUM/2+frame),'.mat']);
    AngMap = zeros(size(f1));
    for index = 1:c
        aaa = f1(:,index);
        fringe =  aaa-s-r;
        y = fringe;
        Hy = hilbert(y);
        HyAng = unwrap(angle(Hy));
        HyAng1 = HyAng - HyAng(T1);
        aaa = f2(:,index);
        fringe = aaa-s-r;
        y = fringe;
        Hy = hilbert(y);
        HyAng = unwrap(angle(Hy));
        HyAng2 = HyAng - HyAng(T1);
        MappingAng = HyAng2 - HyAng1;
        MappingAng = MappingAng/MappingAng(T2);
        AngMap(:,index) = MappingAng;
    end
    MA(:,frame) = mean(AngMap,2);
end
MAmean = mean(MA(:,1:DISPERSION_NUM/2),2);
xi = (linspace(MAmean(1),MAmean(end),pn))';
parfor frame = 1:DISPERSION_NUM
    display(['Calculate CP :',num2str(frame)]);
    f = importdata(['f',num2str(frame),'.mat']);
    CPhaseArray = zeros(size(f));
    for index = 1:c
        aaa = f(:,index);
        fringe = aaa-s-r;
        y = fringe;
        yC = interp1(MAmean,y,xi);
        FyC = abs(fft(yC,NFFT));
        yCAng = unwrap(angle(hilbert(yC)));
        p = polyfit((0:(length(yCAng)-1))',yCAng,1);
        CPhase = yCAng - polyval(p,(0:(length(yCAng)-1))');
        CPhaseArray(:,index) = CPhase;
    end
    CP(:,frame) = mean(CPhaseArray,2);
end
CPhase = mean(CP(:,1:DISPERSION_NUM),2);
CArray = exp(-(1i.*CPhase));
CalStru.MAmean = MAmean;
CalStru.xi = xi;
CalStru.CArray = CArray;
save('imfirst','CalStru');

%% PLOT PSF versus depth
x = linspace(1,NFFT/2,NFFT/2)*x_scale_rate;
for frame=1:PLOT_NUM
    display(['Plot PSF :',num2str(frame)]);
    f=importdata(['f',num2str(frame),'.mat']);
    aaa = f(:,134);
    fringe = aaa-s-r;
    y = fringe;
    yC = interp1(MAmean,y,xi);
    FyC = abs(fft(yC.*CArray,NFFT));
    plot(x,FyC(1:NFFT/2));
    hold on
end
psf_name1 = strcat('psf vs depth','.fig');
psf_name2 = strcat('psf vs depth','.jpg');
saveas(gcf, psf_name1);
saveas(gcf, psf_name2);
end