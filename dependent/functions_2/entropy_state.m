function [ent] = entropy_state(p)
%ENTROPY_STATE Summary of this function goes here
%   Detailed explanation goes here
ent=zeros(size(p,1),1);
for i=1:size(p,1)
    for j=1:size(p,2)
        if p(i,j)~=0
            ent(i,1) = ent(i,1)+(-p(i,j)*log(p(i,j)));
        end
    end
end

end

