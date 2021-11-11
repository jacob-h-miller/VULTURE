%Carbon number-weighted property of a mixture blend. Originally written by 
%Hannah Nguyen for flash point, Jake Miller 3/24/2021
function [Property] =CNumberBlend(CNumber,PropVector,MoleFracs)
wfp = 1./(CNumber.*CNumber);
wfpBlend = wfp.*MoleFracs;
PropBlend= PropVector.*wfpBlend;
Property = sum(PropBlend)/sum(wfpBlend);
end