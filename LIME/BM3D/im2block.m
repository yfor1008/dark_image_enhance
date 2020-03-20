% ͼ��ֿ飬�������任��Ϊ�����ƿ���׼��
% k:���С��p�����ƶ�������lambda2D��delta ������ֵ
% block ���صĿ飬transform_block �任��Ŀ�
% block2row_idx,block2col_idx Ϊ����� ������ �� �����Ͻ���ͼ��������� ��Ӧ��ϵ
function [block,transform_block,block2row_idx,block2col_idx] =im2block(img,k,p,lambda2D,delta)
    [row,col]=size(img);
    %�����ֵ����ʲô��ʽ�أ���
    thres=lambda2D*delta*sqrt(2*log(row*col));
    % r_num:�з��� �� Ӧ���� ���ٸ���
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
            %����ʲô�任�أ���
            tr_b=fft2(block(:,:,cnt));
%             tr_b=dct2(block(:,:,cnt));
            idx=find(abs(tr_b)<thres);
            tr_b(idx)=0;
            transform_block(:,:,cnt)=tr_b;
            cnt=cnt+1;
        end
    end