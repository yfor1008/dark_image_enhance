function [enhanced] = homo_filter(im)
% homo_filter - 同态滤波增强
%
% input:
%   - im: h*w*3, rgb图像
%   - high: float, 高频增益, >1
%   - low: float, 低频增益, [0, 1]
%   - C: float, 锐化系数
%   - sigma: float, 截止频率, 越大图像月亮
% output:
%   - enhanced: h*w*3, 增强后图像
%
% docs:
%   - https://www.cnblogs.com/Imageshop/p/9766056.html
%   - 建议参数: high=2, low=0.2, C=0.1, sigma=max(h, w);
%   - 对正常图像可能存在问题
%

hsv = rgb2hsv(im);
V = round(hsv(:,:,3) * 255);

[h, w] = size(V);
pix_num = h * w;
cx = floor(w/2); % 中心点
cy = floor(h/2);

% 根据图像中亮度低的像素的比例确定是否需要增强
dark_pix_num = length(find(V < 50));
if dark_pix_num < pix_num * 0.05
    enhanced = im;
    return;
end

sigma = max(h, w);
high = 2;
low = 0.2;
C = 0.1;

H = zeros(h, w);
for Y = 1:h
    dist1 = (Y - cy) * (Y - cy);
    for X = 1:w
        dist2 = (X - cx) * (X - cx);
        dist = dist1 + dist2; % 到中心的距离
        H(Y,X) = (high - low) * (1 - exp(-C * (dist / (2 * sigma * sigma)))) + low; % 同态滤波函数
    end
end
H = ifftshift(H); % 反中心化处理

im_log = log(V + 1);
im_fft = fft2(im_log);
im_fft = H .* im_fft; % 滤波
im_ifft = ifft2(im_fft);

V1 = exp(im_ifft) - 1; % 取指数

% 归一化到[0, 1]
val_min = min(V1(:));
val_max = max(V1(:));
range = val_max - val_min;
V1 = (V1 - val_min) / range;
V1 = real(V1);

hsv(:,:,3) = V1;
enhanced = uint8(hsv2rgb(hsv) * 255);

end