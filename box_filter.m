function [filtered] = box_filter(im, win_size)
% box_filter - �����˲�, �൱�ھ�ֵ�˲�
%
% input:
%   - im: h*w, ���˲�ͼ��, ��Ϊgrayͼ��
%   - win_size: int, �˲����ڰ뾶
% output:
%   - filtered: h*w, �˲���ͼ��
%
% docs:
%   - ʹ�û���ͼ���м���
%   - im_dst(y,x)��ʾ��(y,x)���� win_size �뾶�����ڵ��������غ�,
%       ��im_dst(y,x)=sum(im(y-win_size:y+win_size, x-win_size:x+win_size))
%

[h, w, ~] = size(im);
im_dst = zeros(h, w);

% ���м���
im_row = cumsum(im, 1);
im_dst(1:win_size+1, :) = im_row(1+win_size:win_size*2+1, :); % ��Եwin_size+1��/��
im_dst(win_size+2:h-win_size, :) = im_row(win_size*2+2:h, :) - im_row(1:h-win_size*2-1, :); % �м�����
im_dst(h-win_size+1:h, :) = repmat(im_row(h, :), [win_size, 1]) - im_row(h-win_size*2:h-win_size-1, :); % ��Ե����

% ���м���
im_col = cumsum(im_dst, 2);
im_dst(:, 1:win_size+1) = im_col(:, 1+win_size:win_size*2+1);
im_dst(:, win_size+2:w-win_size) = im_col(:, win_size*2+2:w) - im_col(:, 1:w-win_size*2-1);
im_dst(:, w-win_size+1:w) = repmat(im_col(:, w), [1, win_size]) - im_col(:, w-win_size*2:w-win_size-1);

filtered = im_dst;

end