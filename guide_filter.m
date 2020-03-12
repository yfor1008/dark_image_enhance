function [filtered] = guide_filter(im, guide_im, win_size, eps)
% guide_filter - �����˲�
%
% input:
%   - im: h*w, ���˲�ͼ��, ��Ϊgrayͼ��
%   - guide_im: h*w*c, ����ͼ��, c=3ʱΪrgbͼ��, c=1ʱΪgrayͼ��
%   - win_size: int, �˲����ڰ뾶
%   - eps: float, Խ��ͼ��Խģ��
% output:
%   - filtered: h*w, �˲���ͼ��
%
% docs:
%

win_size = win_size * 2 + 1;

im = double(im);
guide_im = double(guide_im);

[h, w, c] = size(guide_im);

N = box_filter(ones(h, w), win_size); % �˲����������ظ���

if c == 1
    mean_I = box_filter(im, win_size) ./ N;
    mean_G = box_filter(guide_im, win_size) ./ N;
    mean_IG = box_filter(im .* guide_im, win_size) ./ N;
    cov_IG = mean_IG - mean_I .* mean_G; % Э����
    
    mean_GG = box_filter(guide_im .* guide_im, win_size) ./ N;
    var_G = mean_GG - mean_G .* mean_G;

    a = cov_IG ./ (var_G + eps); % ��ʽ5
    b = mean_I - a .* mean_G; % ��ʽ6

    mean_a = box_filter(a, win_size) ./ N;
    mean_b = box_filter(b, win_size) ./ N;

    filtered = mean_a .* guide_im + mean_b; % ��ʽ8

elseif c == 3
    mean_Gr = box_filter(guide_im(:,:,1), win_size) ./ N;
    mean_Gg = box_filter(guide_im(:,:,2), win_size) ./ N;
    mean_Gb = box_filter(guide_im(:,:,3), win_size) ./ N;

    mean_I = box_filter(im, win_size) ./ N;

    mean_Gr_I = box_filter(guide_im(:,:,1) .* im, win_size) ./ N;
    mean_Gg_I = box_filter(guide_im(:,:,2) .* im, win_size) ./ N;
    mean_Gb_I = box_filter(guide_im(:,:,3) .* im, win_size) ./ N;

    cov_Gr_I = mean_Gr_I - mean_Gr .* mean_I;
    cov_Gg_I = mean_Gg_I - mean_Gg .* mean_I;
    cov_Gb_I = mean_Gb_I - mean_Gb .* mean_I;

    var_Grr_I = box_filter(guide_im(:,:,1) .* guide_im(:,:,1), win_size) ./ N - mean_Gr .* mean_Gr;
    var_Grg_I = box_filter(guide_im(:,:,1) .* guide_im(:,:,2), win_size) ./ N - mean_Gr .* mean_Gg;
    var_Grb_I = box_filter(guide_im(:,:,1) .* guide_im(:,:,3), win_size) ./ N - mean_Gr .* mean_Gb;
    var_Ggg_I = box_filter(guide_im(:,:,2) .* guide_im(:,:,2), win_size) ./ N - mean_Gg .* mean_Gg;
    var_Ggb_I = box_filter(guide_im(:,:,2) .* guide_im(:,:,3), win_size) ./ N - mean_Gg .* mean_Gb;
    var_Gbb_I = box_filter(guide_im(:,:,3) .* guide_im(:,:,3), win_size) ./ N - mean_Gb .* mean_Gb;

    a = zeros(h, w, 3);
    for row = 1:h
        for col = 1:w
            sigma = [var_Grr_I(row, col), var_Grg_I(row, col), var_Grb_I(row, col);
                     var_Grg_I(row, col), var_Ggg_I(row, col), var_Ggb_I(row, col);
                     var_Grb_I(row, col), var_Ggb_I(row, col), var_Gbb_I(row, col)];
            sigma = sigma + eps * eye(3);

            cov_IG = [cov_Gr_I(row, col), cov_Gg_I(row, col), cov_Gb_I(row, col)];

            a(row, col, :) = cov_IG * inv(sigma + eps * eye(3));
        end
    end

    b = mean_I - a(:,:,1) .* mean_Gr - a(:,:,2) .* mean_Gg - a(:,:,3) .* mean_Gb;

    filtered = (box_filter(a(:,:,1), win_size) .* guide_im(:,:,1) ...
              + box_filter(a(:,:,2), win_size) .* guide_im(:,:,2) ...
              + box_filter(a(:,:,3), win_size) .* guide_im(:,:,3) ...
              + box_filter(b, win_size)) ./ N;
    
else
    error('channel of im must be 1 or 3 !');
end

end