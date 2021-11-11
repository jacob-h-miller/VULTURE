% Determine the boiling point of a mixture of species using Raoult's Law,
% assuming ideal gas and liquid phases.
%  Jake Miller, June 18, 2021
function Tboil = RaoultsLawCalc(molfracs,VPsubset)
VPs = zeros(size(VPsubset,1),1);
i_save = 0;
molfracsnew = molfracs ./ sum(molfracs);
for i=1:size(VPs,1)
    VPs(i) = dot(molfracsnew,VPsubset(i,:));
    if VPs(i) > 1 % stop calculating if we've hit vapor pressure=1 atm (boiling)
        i_save = i;
        break
    end
end

T_low = -43; % Temp, in C, of lowest vapor pressure data.
Tboil = T_low + i_save; % Boiling temperature, in C.

end
