% Check to see if a fuel meets all specified criteria
function yesorno = checkspecs(FuelSpecs, SpecLimits)
yesorno = zeros(size(FuelSpecs));
for i=1:size(FuelSpecs,2)
    if (FuelSpecs(i) > SpecLimits(1,i)) && (FuelSpecs(i) < SpecLimits(2,i))
        yesorno(i) = 1; % Fuel meets that specific spec
    end
end