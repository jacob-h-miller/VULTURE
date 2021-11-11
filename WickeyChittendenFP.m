% Calculate the flash point of a liquid mixture using the Wickey and
% Chittenden model

function FP = WickeyChittendenFP(volfracs,BlendIndices)

BlendIndexmix = dot(volfracs, BlendIndices);
FP_Kelvin = 2414/(6.1188+log10(BlendIndexmix))-230.56; % Flash point in Kelvin
FP = FP_Kelvin - 273; % Flash point in celsius

end