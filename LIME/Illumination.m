function T = Illumination(L, method)
% Inputs:
%        L:
% Outputs:
%        T:

if strcmp(method, 'max_c') == 1
    T = max(L, [], 3);
elseif strcmp(method, 'min_c') == 1
    T = min(L, [], 3); 
end 

end