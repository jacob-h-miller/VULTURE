% Figure out if species to cut to meet spec comes from the top or bottom
% Jake Miller, April 14, 2021
function cut_direction = Top_Bottom_Cut(Spec_Failed, Low_Hi)
if (Spec_Failed == "Boiling Point") || (Spec_Failed == "Flash Point") || (Spec_Failed == "Lower Heating Value") ||...
        (Spec_Failed == "Viscosity") || (Spec_Failed == "Melting Point") || (Spec_Failed == "Cetane Number")
    if Low_Hi == "High"
        cut_direction = "Bottom";
    elseif Low_Hi == "Low"
        cut_direction = "Top";
    else cut_direction = "Low_Hi Error";
    end
elseif (Spec_Failed == "Water Solubility") || (Spec_Failed == "Octane Number")
    if Low_Hi == "High"
        cut_direction = "Top";
    elseif Low_Hi == "Low"
        cut_direction = "Bottom";
    else cut_direction = "Low_Hi Error";
    end
else cut_direction = "Failed Spec Definition Error"
end
end