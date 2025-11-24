
function [tsw] = transition_subject_wise(posterior_states_evolution,list)
% posterior_states_evolution is state evolution map for each subject,...
%  ...defined as a 1*tr vectors extracted from estStatesCell from original code.
% list is detected states during group-level BSDS, sorted by their oocupancy.
% [tsw] stored all results as a matrix.
% basicly this code employed string-matching stragety to count through whole session. 
tsw=zeros(0);
 for i=1:length(list)   % counting by row, e.g. s1 to s1, s1 to s2,...
  for k=1:length(list)
  b=[num2str(list(i)),num2str(list(k))];   % matching string
  count=0;
   for j=1:length(posterior_states_evolution)-1
     if strcmp([strcat(num2str(posterior_states_evolution(j)),num2str(posterior_states_evolution(j+1)))],b)
     count=count+1;
     end
   end
  tsw(i,k)=count;
  end
   if find(posterior_states_evolution()==list(i))~=0;  % make sure denominator is not equal to 0;
       tsw(i,:)= tsw(i,:)/sum(tsw(i,:));
    end
  end
 end