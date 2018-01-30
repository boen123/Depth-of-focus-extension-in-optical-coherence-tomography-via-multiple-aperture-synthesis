function [] = print_image(name,img)
% this code aims to plot images

img_temp = abs(img);
D1_step1 = (img_temp-min(min(img_temp)))./(max(max(img_temp))-min(min(img_temp)))*255;
dir_set(4);
fname1 = sprintf(name);
imwrite(uint8(D1_step1),fname1,'BMP');
end