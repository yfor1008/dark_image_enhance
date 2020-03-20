% 3DÄæ±ä»»
function [blk_est]=inv_transform_3d(blk_tran3d,tran_mode)
    global blk_tran1d_s;
    global blk_2d_s;
    [m,n,blk_num]=size(blk_tran3d);
 
    blk_invtran1d=zeros(m,n,blk_num);
    blk_est=zeros(m,n,blk_num);
 
    if(tran_mode==0)    %fft
        for i=1:1:m
            for j=1:1:n
                blk_invtran1d(i,j,:)=ifft(blk_tran3d(i,j,:));
            end
        end
        for i=1:1:blk_num
            blk_est(:,:,i)=ifft2(blk_invtran1d(:,:,i));
        end
    elseif(tran_mode==1)  %dct
        for i=1:1:m
            for j=1:1:n
                blk_invtran1d(i,j,:)=idct(blk_tran3d(i,j,:));
            end
        end
        for i=1:1:blk_num
            blk_est(:,:,i)=idct2(blk_invtran1d(:,:,i));
        end
    elseif(tran_mode==2)    %dwt
        blk_num=length(blk_2d_s);
        blk_c=waverec2(blk_tran3d,blk_tran1d_s,'haar');
        blk_est=[];
        for i=1:1:blk_num
            blk_est(:,:,i)=waverec2(blk_c(:,i),blk_2d_s{i},'Bior1.5');
        end
 
    else
        error('tran_mode error');
    end