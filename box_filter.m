function [filtered] = box_filter(im, win_size)
% box_filter - 盒子滤波, 相当于均值滤波
%
% input:
%   - im: h*w, 待滤波图像, 需为gray图像
%   - win_size: int, 滤波窗口半径
% output:
%   - filtered: h*w, 滤波后图像
%
% docs:
%   - 使用积分图进行加速
%   - im_dst(y,x)表示点(y,x)在以 win_size 半径邻域内的所有像素和,
%       即im_dst(y,x)=sum(im(y-win_size:y+win_size, x-win_size:x+win_size))
%

[h, w, ~] = size(im);
im_dst = zeros(h, w);

% 按行计算
im_row = cumsum(im, 1);
im_dst(1:win_size+1, :) = im_row(1+win_size:win_size*2+1, :); % 边缘win_size+1行/列
im_dst(win_size+2:h-win_size, :) = im_row(win_size*2+2:h, :) - im_row(1:h-win_size*2-1, :); % 中间数据
im_dst(h-win_size+1:h, :) = repmat(im_row(h, :), [win_size, 1]) - im_row(h-win_size*2:h-win_size-1, :); % 边缘数据

% 按列计算
im_col = cumsum(im_dst, 2);
im_dst(:, 1:win_size+1) = im_col(:, 1+win_size:win_size*2+1);
im_dst(:, win_size+2:w-win_size) = im_col(:, win_size*2+2:w) - im_col(:, 1:w-win_size*2-1);
im_dst(:, w-win_size+1:w) = repmat(im_col(:, w), [1, win_size]) - im_col(:, w-win_size*2:w-win_size-1);

filtered = im_dst;

end