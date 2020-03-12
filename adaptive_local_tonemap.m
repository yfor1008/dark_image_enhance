function [enhanced] = adaptive_local_tonemap(im)
% adaptive_local_tonemap - ����Ӧ�ֲ�����
%
% input:
%   - im: h*w*3, rgbͼ��
% output:
%   - enhanced: h*w*3, ��ǿ��ͼ��
%
% docs��
%   - Adaptive Logarithmic Mapping For Displaying High Contrast Scenes
%   - https://www.cnblogs.com/Imageshop/p/8420828.html
%   - https://blog.csdn.net/fb_941219/article/details/88900956
%   - https://blog.csdn.net/just_sort/article/details/84066390
%   - opencv��ʵ��: https://docs.opencv.org/3.0-beta/modules/photo/doc/hdr_imaging.html
%

[h,w,~] = size(im);
im = double(im);
enhanced = zeros(h, w, 3);

im = im / 255;

% RGBתXYZ
X = (0.4124 * im(:,:,1) + 0.3576 * im(:,:,2) + 0.1805 * im(:,:,3));
Y = (0.2126 * im(:,:,1) + 0.7152 * im(:,:,2) + 0.0722 * im(:,:,3));
Z = (0.0193 * im(:,:,1) + 0.1192 * im(:,:,2) + 0.9505 * im(:,:,3));

bias = 0.75;
Lwmax = max(Y(:));
Lwav = mean(Y(:));
Lwmax = Lwmax / Lwav;
bias_p = log(bias)/log(0.5); % ��ʽ3�е�ָ��
div = log10(Lwmax+1); % ��ʽ4�еĳ���

Y_new = log(Y/Lwav+1) / div ./ log(2+power(Y/Lwav/Lwmax, bias_p)*8); % ��ʽ4
X_new = Y_new ./ Y .* X; % ͬ��������X/Zͨ��
Z_new = Y_new ./ Y .* Z;

% XYZתRGB
enhanced(:,:,1) =  3.2410 * X_new - 1.5374 * Y_new - 0.4986 * Z_new;
enhanced(:,:,2) = -0.9692 * X_new + 1.8760 * Y_new + 0.0416 * Z_new;
enhanced(:,:,3) =  0.0556 * X_new - 0.2040 * Y_new + 1.0570 * Z_new;

% % gamma correction
% enhanced = gamma_correction(enhanced);

enhanced = uint8(enhanced*255);

end

function [corrected] = gamma_correction(x)
% gamma_correction - gamma����
%

corrected = x;
idx1 = (corrected <= 0.05);
idx2 = (corrected > 0.05);
corrected(idx1) = corrected(idx1) * 2.64;
corrected(idx2) = 1.099 * power(corrected(idx2), 0.9/2.2) - 0.099;

end