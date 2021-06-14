function [y,idx] = nanmin(a,b,dim)


if nargin==1
    [y,idx] = min(a);
end

if nargin==2
    [y,idx] = min(a,b);
end

if nargin==3
    if isempty(b)
        [y,idx] = min(a,[],dim); 
    else
        error('"b" have to be empty')
    end
end