
function ab = mean_edge(ginv,Hinv)
NElement=max(size(Hinv));

ed_inv= zeros(6*NElement,1);  

for ii=1:NElement

   elind = Hinv(ii,:);

   gg = ginv(elind,:);

   %edges

   temp = [norm(gg(2,:)-gg(1,:)); norm(gg(3,:)-gg(1,:));norm(gg(4,:)-gg(1,:)); norm(gg(3,:)-gg(2,:)); norm(gg(4,:)- gg(2,:)); norm(gg(4,:)-gg(3,:))];
 
   ed_inv(6*(ii-1)+1:6*ii) = temp;
   ab(ii) = mean(temp);
 

%    if(max(temp)>=3)
% 
%        keyboard
% 
%    end

end 

ab = min(ab);