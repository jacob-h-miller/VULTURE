% Converts a vector of (alcohol) mass fractions (massfracs) to a corresponding vector
% of (alcohol) volume fractions (volfrac) using the species densities (rho).
function volfracs = volfraccalc(massfracs, rho)
vtot = massfracs .* rho;
volcalc = massfracs ./ rho;
volfracs = volcalc ./ sum(volcalc);
end