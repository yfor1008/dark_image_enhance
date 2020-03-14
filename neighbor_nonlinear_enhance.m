function [enhanced] = neighbor_nonlinear_enhance(im)
% neighbor_nonlinear_enhance - 邻域非线性增强
%
% input:
%   - im: h*w*3, rgb图像
% output:
%   - enhanced: h*w*3, 增强后图像
%
% docs:
%   - 《An Integrated Neighborhood Dependent Approach for Nonlinear Enhancement of Color Images》
%   - https://blog.csdn.net/weixin_44690935/article/details/102680513
%

im = im2double(im); % [0, 1]
I1 = rgb2gray(im);

enhanced = zeros(size(im));

In = (I1 .^ 0.24 + (1 - I1) * 0.5 + I1 .^ 2) / 2;

sigma1 = 5;
window1 = double(uint8(sigma1 * 3) * 2 + 1);
sigma2 = 20;
window2 = double(uint8(sigma2 * 3) * 2 + 1);
sigma3 = 240;
window3 = double(uint8(sigma3 * 3) * 2 + 1);

im_g1 = imgaussfilt(I1, sigma1, 'FilterSize', window1);
im_g2 = imgaussfilt(I1, sigma2, 'FilterSize', window2);
im_g3 = imgaussfilt(I1, sigma3, 'FilterSize', window3);

r1 = im_g1 ./ I1;
R1 = In .^ r1;

r2 = im_g2 ./ I1;
R2 = In .^ r2;

r3 = im_g3 ./ I1;
R3 = In .^ r3;

R = (R1 + R2 + R3)/3;

sc = min(R ./ I1, 4); % 增加限制
enhanced(:,:,1) = sc .* im(:,:,1);
enhanced(:,:,2) = sc .* im(:,:,2);
enhanced(:,:,3) = sc .* im(:,:,3);

enhanced = uint8(enhanced * 255);

end