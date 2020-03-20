% �������ƿ�
% ik,jk: ������Ͻ� ���ꡣ
%bn_r: block_num_row,�ֿ�֮��,�����еķ������ж��ٸ��顣
%�� tran_block ������룬�õ���sim_blk ȴ������block.
%np: �������ڵĴ�С/���ƶ����ࡣ
% tau: ������ƾ���
% max_sim_num: ������ƿ� ����
function [sim_blk,sim_num,sim_blk_idx]=search_similar_block(ik,jk,block,tran_block,k,np,bn_r,bn_c,tau,max_sim_num)
    % �������ڵ� ���Ͻ� �� ���½� �Ŀ�����ꡣs,e:start,end.
    in_s=max(ik-floor(np/2),1);
    jn_s=max(jk-floor(np/2),1);
    in_e=min(ik+floor(np/2),bn_r);
    jn_e=min(jk+floor(np/2),bn_c);
    % ��ǰ�ο���
    ref_blk=tran_block(:,:,((ik-1)*bn_c+jk));
    k2=k*k;
    cnt=0;
     
    %dist=[];
    %blk_idx=[];
    %����ο�����ͼ��ı�Ե��������Χ����������Ͳ�����������n*n
    %���Բ��ܵ�������cnt ���Ƴ�idx.
    %for ii=in_s:1:in_e
    %   for jj=jn_s:1:jn_e
    %       cnt=cnt+1;
    %       %idx=ii*bn_c+jj;
    %       idx=(ii-1)*bn_c+(jj-1)+1;
    %       cur_blk=tran_block(:,:,idx);
    %       blk_idx(cnt)=idx;
    %       dist(cnt)=norm(cur_blk-ref_blk);% ֻ�ǱȽϴ�С����û��Ҫ��һ���� /k2;
    %   end
    %end
    %���������ƿ�ķ���Ҫ������Ŀ�һЩ
 
    ii=in_s:1:in_e;
    jj=jn_s:1:jn_e;
    [II,JJ]=meshgrid(ii,jj);
    IDX=(II-1)*bn_c+JJ;
    blk_idx=IDX(:);
    cur_blk=tran_block(:,:,blk_idx);
    cnt=size(cur_blk,3);
    ref_blk_mat=repmat(ref_blk,[1,1,cnt]);
    % shit! Ҫ��Ҫ��norm ? ����ֵ�ܺ�����������Ƴ̶ȣ�
    % norm(cur_blk-ref_blk_mat);
    delta_blk=cur_blk-ref_blk_mat;
    dist=sum(sum(delta_blk.*delta_blk,1),2);
 
    % �����ҵ������ƵĿ�
    [dist_sort,dist_idx]=sort(dist);
    max_num=min(cnt,max_sim_num);%���ܵ�ǰ���и�����û��max_sim_num��
    if(dist_sort(max_num)<tau)
        sim_num=max_num;
    else
        sim_num=sum(dist_sort(1:max_num)<tau);
    end
    cnt_idx=dist_idx(1:sim_num);
    sim_blk_idx=blk_idx(cnt_idx);
    sim_blk=block(:,:,sim_blk_idx);