% Converts a vector of (alcohol) mole fractions (molfracs) to a corresponding vector
% of (alcohol) weight fractions (wtfrac) using the species molecular
% weights (MWs).
function wtfrac = moltowtfrac(molfracs,MWs)
masses = molfracs .* MWs;
wtfrac = masses ./ sum(masses);
end