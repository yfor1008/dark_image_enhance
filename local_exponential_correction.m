function [enhanced] = local_exponential_correction(im)
% local_exponential_correction - 局部指数矫正
%
% input:
%   - im: h*w*3, rgb图像
% output:
%   - enhanced: h*w*3, 增强后图像
%
% docs:
%   - 《Contrast image correction method》
%   - 
%

luminance = rgb2gray(im);
[h, w] = size(luminance);
% pix_num = h * w;

im = double(im);
enhanced = zeros(h, w, 3);

mask = guide_filter(luminance, luminance, max(round(max(h, w)*0.01), 5), 25);
% mask = imgaussfilt(luminance, 'FilterSize', max(round(max(h, w)*0.01), 5));
luminance = double(luminance);
mask = 255 - mask;
mask = double(mask);
% luminance = double(luminance);

new_lum = power(luminance/255, power(2, (128-mask)/128.0));
new_lum = round(new_lum * 255);

for Y = 1:h
    for X = 1:w
        old = luminance(Y,X);
        if old == 0
            enhanced(Y,X,:) = [0,0,0]';
        else
            new = new_lum(Y,X);
            enhanced(Y,X,:) = (new * im(Y,X,:) / old + im(Y,X,:) + new - old) / 2;
        end
    end
end
enhanced = uint8(enhanced);

end