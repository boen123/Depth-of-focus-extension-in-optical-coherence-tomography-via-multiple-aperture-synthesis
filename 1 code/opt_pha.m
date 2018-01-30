function [opt_pha1, opt_pha2, opt_pha3, opt_pha4, opt_pha5] = opt_pha(img1, img2, img3, img4, img5, pha_num)
% Achieve the optimal phase to correct the defocus-induced effect using the fourier shift theorem.

phi = linspace(0,2*pi,pha_num+1);
ASS = 0; 
for num1 = 1:pha_num
    im1 = img1 * exp(-1i * phi(num1));
    for num2 = 1:pha_num
        im2 = img2 * exp(-1i * phi(num2));
        for num3 = 1:pha_num
            im3 = img3 * exp(-1i * phi(num3));
            for num4 = 1:pha_num
                im4 = img4 * exp(-1i * phi(num4));
                for num5 = 1:pha_num
                    im5 = img5 * exp(-1i * phi(num5));
                    %% calculate the max value of the coherent sum (complex sum). Optimal corrected phase obtained when ASS is maximum.
                    sum_img = abs(im1 + im2 + im3 + im4 + im5);
                    ASS_tmp = max(max(sum_img));
                    if ASS_tmp >= ASS
                        ASS = ASS_tmp;
                        opt_pha1 = phi(num1);
                        opt_pha2 = phi(num2);
                        opt_pha3 = phi(num3);
                        opt_pha4 = phi(num4);
                        opt_pha5 = phi(num5);
                    end
                end
            end
        end
    end
end
