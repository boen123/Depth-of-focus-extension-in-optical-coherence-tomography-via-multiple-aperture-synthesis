function img_sfted = fsft(in_img, opt_axi_sft, opt_lat_sft, CArray, NFFT)
% Shifts input image using the fourier shift theorem.

[M, N] = size(in_img);
z_sft = exp(+1i * 2 * pi * opt_axi_sft * [0:floor(M/2)-1 floor(-M/2):-1]' / M);
x_sft = exp(+1i * 2 * pi * opt_lat_sft * [0:floor(N/2)-1 floor(-N/2):-1] / N);   
img_sfted = fft((in_img .*CArray) .* (z_sft * x_sft),NFFT);   
end