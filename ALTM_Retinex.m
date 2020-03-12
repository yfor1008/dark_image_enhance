function outval = ALTM_Retinex(I)
% ALTM_Retinex - ����Ӧ�ֲ�����
%
% input:
%   - I: h*w*3, rgbͼ��
% output:
%   - outval: h*w*3, ��ǿ��ͼ��
%
% docs:
%   - ��Adaptive Local Tone Mapping Based on Retinex for High Dynamic Range Images��
%   - https://www.cnblogs.com/Imageshop/p/9460334.html
%

II = im2double(I);
Ir = double(II(:,:,1));
Ig = double(II(:,:,2));
Ib = double(II(:,:,3));

% Global Adaptation
Lw = 0.299 * Ir + 0.587 * Ig + 0.114 * Ib; % input world luminance values
Lwmax = max(max(Lw)); % the maximum luminance value
[m, n] = size(Lw);
Lwaver = exp(sum(sum(log(0.001 + Lw))) / (m * n)); % ��ʽ5, log-average luminance
Lg = log(Lw / Lwaver + 1) / log(Lwmax / Lwaver + 1); % ��ʽ4
gain = Lg ./ Lw;
gain(Lw == 0) = 0;
outval = cat(3, gain .* Ir, gain .* Ig, gain .* Ib);

end