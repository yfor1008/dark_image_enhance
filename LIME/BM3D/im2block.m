% 图像分块，并且做变换，为找相似块做准备
% k:块大小，p：块移动步长，lambda2D，delta 收缩阈值
% block 返回的块，transform_block 变换后的块
% block2row_idx,block2col_idx 为保存的 块索引 与 块左上角在图像中坐标的 对应关系
function [block,transform_block,block2row_idx,block2col_idx] =im2block(img,k,p,lambda2D,delta)
    [row,col]=size(img);
    %这个阈值该用什么公式呢？？
    thres=lambda2D*delta*sqrt(2*log(row*col));
    % r_num:行方向 上 应该有 多少个块
    r_num=floor((row-k)/p)+1;
    c_num=floor((col-k)/p)+1;
    block=zeros(k,k,r_num*c_num);
    block2row_idx=[];
    block2col_idx=[];
    cnt=1;
    for i=0:1:r_num-1
        rs=1+i*p;
        for j=0:1:c_num-1
            cs=1+j*p;
            block(:,:,cnt)=img(rs:rs+k-1,cs:cs+k-1);
            block2row_idx(cnt)=rs;
            block2col_idx(cnt)=cs;
            %该用什么变换呢？？
            tr_b=fft2(block(:,:,cnt));
%             tr_b=dct2(block(:,:,cnt));
            idx=find(abs(tr_b)<thres);
            tr_b(idx)=0;
            transform_block(:,:,cnt)=tr_b;
            cnt=cnt+1;
        end
    end