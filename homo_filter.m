function [enhanced] = homo_filter(im)
% homo_filter - ̬ͬ�˲���ǿ
%
% input:
%   - im: h*w*3, rgbͼ��
%   - high: float, ��Ƶ����, >1
%   - low: float, ��Ƶ����, [0, 1]
%   - C: float, ��ϵ��
%   - sigma: float, ��ֹƵ��, Խ��ͼ������
% output:
%   - enhanced: h*w*3, ��ǿ��ͼ��
%
% docs:
%   - https://www.cnblogs.com/Imageshop/p/9766056.html
%   - �������: high=2, low=0.2, C=0.1, sigma=max(h, w);
%   - ������ͼ����ܴ�������
%

hsv = rgb2hsv(im);
V = round(hsv(:,:,3) * 255);

[h, w] = size(V);
pix_num = h * w;
cx = floor(w/2); % ���ĵ�
cy = floor(h/2);

% ����ͼ�������ȵ͵����صı���ȷ���Ƿ���Ҫ��ǿ
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
        dist = dist1 + dist2; % �����ĵľ���
        H(Y,X) = (high - low) * (1 - exp(-C * (dist / (2 * sigma * sigma)))) + low; % ̬ͬ�˲�����
    end
end
H = ifftshift(H); % �����Ļ�����

im_log = log(V + 1);
im_fft = fft2(im_log);
im_fft = H .* im_fft; % �˲�
im_ifft = ifft2(im_fft);

V1 = exp(im_ifft) - 1; % ȡָ��

% ��һ����[0, 1]
val_min = min(V1(:));
val_max = max(V1(:));
range = val_max - val_min;
V1 = (V1 - val_min) / range;
V1 = real(V1);

hsv(:,:,3) = V1;
enhanced = uint8(hsv2rgb(hsv) * 255);

end