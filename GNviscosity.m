% Viscosity of liquid mixtures using the Grunberg-Nissan Equation with no
% interaction parameter: ln(eta_mix)=sum(x_i*ln(eta_i)), where etas are
% viscosities and x_is are mole fractions. Jake Miller, 6/18/2021.
% Viscosities are delivered and reported in cSt. Source of the equation: https://link.springer.com/content/pdf/10.1007%2Fs10765-011-1100-1.pdf

function mu = GNviscosity(mole_fracs, viscosities)

logsum = 0;
for i=1:size(mole_fracs,1)
    contribution_i = mole_fracs(i) * log(viscosities(i));
    if viscosities(i) == 0
        contribution_i = 0;
    end
    logsum = logsum + contribution_i;
end

mu = exp(logsum);

if isnan(mu)
    disp('NaN viscosities');
end

end