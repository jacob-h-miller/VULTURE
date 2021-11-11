%Linear blending for a given product (assumption). Jake Miller, 3/24/2021.
function [Property]=LinearBlend(PropVector,WtFractions)
PropBlend = PropVector.*WtFractions;
Property = sum(PropBlend)/sum(WtFractions);
end
