function dy = ODEs_ketone (time,y);
global k n total 

% y represents all components
% n is the number of acid reactants
% n(n+1)/2 is the number of ketone product and the number of ketonization
% rxn

dy=zeros(total,1);
% k is rate constant, 15 reaction
% 1 - C4-C4      % 6 - C5-C5    % 10 - C6-C6   % 13 - C7-C7  % 15 - C8-C8
% 2 - C4-C5      % 7 - C5-C6    % 11 - C6-C7   % 14 - C7-C8
% 3 - C4-C6      % 8 - C5-C7    % 12 - C6-C8    
% 4 - C4-C7      % 9 - C5-C8
% 5 - C4-C8
  
% rates of reaction
% ri 

r = zeros(n,n);
for i = 1:n
    for j = i:n        
        r(j,i)= -k(j,i)*y(i)*y(j);
    end
end
% Component concentration profile
% dy(i)
dy= zeros(total,1);
for i = 1:n
    dy(i)= sum(r(1:n,i)) + sum(r(i,1:n));  
    m = n*(i-1)-i*(i-1)/2;
    for j = i:n
        dy(j+m+n) = -r(j,i);
    end  

end
 