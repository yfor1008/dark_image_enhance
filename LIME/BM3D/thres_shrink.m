 % ��ֵ������absС��thres�� ��ֵ 0
function [val]=thres_shrink(data,thres)
    val=data;
    idx=find(abs(data)<thres);
    val(idx)=0;
end