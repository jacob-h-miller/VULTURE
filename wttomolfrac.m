% Converts a vector of weight fractions (wtfracs) to a corresponding vector
% of mole fractions (molefrac) using the species molecular
% weights (MWs).
function molfracs = wttomolfrac(wtfracs,MWs)
moles = wtfracs ./ MWs;
molfracs = moles ./ sum(moles);
end