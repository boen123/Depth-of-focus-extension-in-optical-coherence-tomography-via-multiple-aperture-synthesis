function [z_opt_sft1, z_opt_sft2, z_opt_sft3, z_opt_sft4, z_opt_sft5] = opt_z_sft(in_img1, in_img2, in_img3, in_img4, in_img5, z_sft_st, z_sft_end, z_st, z_end, x_st, x_end)
% calculate the optimal axial shift using the fourier shift theorem.

[N, M] = size(in_img3);
ASS = 0;
kk = [0:floor(N/2)-1 floor(-N/2):-1]' / N;
tt = [0:floor(M/2)-1 floor(-M/2):-1] / M;
for z_sft_tmp1 = z_sft_st:z_sft_end
    % shift phase term
    sft1_z = exp(+1i * 2 * pi * z_sft_tmp1 * kk);
    sft1_x = exp(+1i * 2 * pi * 0 * tt);
    % shift image
    frq1 = fft(in_img1 .* (sft1_z * sft1_x),4096);
    for z_sft_tmp2 = z_sft_st:z_sft_end
        z_sft_tmp2;
        % shift phase term
        sft2_z = exp(+1i * 2 * pi * z_sft_tmp2 * kk);
        sft2_x = exp(+1i * 2 * pi * 0 * tt);
        % shift image
        frq2 = fft(in_img2 .* (sft2_z * sft2_x),4096);
        for z_sft_tmp3 = z_sft_st:z_sft_end
            z_sft_tmp3;
            sft3_z = exp(+1i * 2 * pi * z_sft_tmp3 * kk);
            sft3_x = exp(+1i * 2 * pi * 0 * tt);
            frq3 = fft(in_img3 .* (sft3_z * sft3_x),4096);
            for z_sft_tmp4 = z_sft_st:z_sft_end
                sft4_z = exp(+1i * 2 * pi * z_sft_tmp4 * kk);
                sft4_x = exp(+1i * 2 * pi * 0 * tt);
                frq4 = fft(in_img4 .* (sft4_z * sft4_x),4096);
                for z_sft_tmp5 = z_sft_st:z_sft_end
                    sft5_z = exp(+1i * 2 * pi * z_sft_tmp5 * kk);
                    sft5_x = exp(+1i * 2 * pi * 0 * tt);
                    frq5 = fft(in_img5 .* (sft5_z * sft5_x),4096);
                    % sum 5 images after axial shift
                    frq_ass = abs(frq1(z_st:z_end,x_st:x_end) + frq2(z_st:z_end,x_st:x_end)...
                        + frq3(z_st:z_end,x_st:x_end) + frq4(z_st:z_end,x_st:x_end)...
                        + frq5(z_st:z_end,x_st:x_end));
                    ASS_temp = max(max(frq_ass));
                    if ASS_temp >= ASS
                        ASS = ASS_temp;
                        z_opt_sft1 = z_sft_tmp1;
                        z_opt_sft2 = z_sft_tmp2;
                        z_opt_sft3 = z_sft_tmp3;
                        z_opt_sft4 = z_sft_tmp4;
                        z_opt_sft5 = z_sft_tmp5;
                    end
                end
            end
        end
    end
end
end
