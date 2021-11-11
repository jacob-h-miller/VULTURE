% Calculate the flash point of a liquid mixture using the Torabian and
% Sobati model #1 (https://doi.org/10.1016/j.psep.2017.07.020)

function FP = TorabianSobatiFP(molfracs,BlendIndices)

BlendIndexmix = dot(molfracs, BlendIndices);
FP_Kelvin = BlendIndexmix ^-0.0246; % Flash point in Kelvin
FP = FP_Kelvin - 273; % Flash point in celsius

end