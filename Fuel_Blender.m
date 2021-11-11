% Mix two fuels until they fall out of a specified specification.
%  Jake Miller, April 14, 2021
function [max_blend, mix_props, FailedProp, HiLo] = Fuel_Blender(bio_props, petro_props, Spec_Limits, blendfrac, blendstock_type, fuel_class, biomatrix)
nspec = size(Spec_Limits,2);
blend_meetspec = 1;
mix_props = [blendfrac, petro_props]; % zero blend fraction is just pure petrofuel!
prop_list = ["Boiling Point","Flash Point","Lower Heating Value","Viscosity","Melting Point","Water Solubility","Cetane Number","Octane Number"];
while blend_meetspec == 1
    blendfrac = blendfrac + 0.01;
    mixspecs = blendmix(blendfrac, bio_props, petro_props, blendstock_type, fuel_class, biomatrix);
    mix_yesno = checkspecs(mixspecs,Spec_Limits);
    if sum(mix_yesno) == nspec
        mix_props = [mix_props; blendfrac, mixspecs];
        if blendfrac == 1
            blend_meetspec = 0; % This blend will meet specs, but it's not a blend anymore! Just exit and return its specs.
        end
    else blend_meetspec = 0;
        max_blend = blendfrac - 0.01; % maximum diesel blend limit is 1% less than we just tried
        mix_props = [mix_props; blendfrac, mixspecs];
        for i=1:nspec
            if mix_yesno(i) == 0
                FailedProp = prop_list(i);
                if mixspecs(i) < Spec_Limits(1,i)
                    HiLo = "Low";
                elseif mixspecs(i) > Spec_Limits(2,i)
                    HiLo = "High";
                else
                    HiLo = "Error";
                end
            end
        end
    end
end
end