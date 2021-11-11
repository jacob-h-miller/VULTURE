%  Figure out the best destinations for secondary blending with "waste"
%  products from initial cuts of a blended alcohol-based biofuel.
%  Jake Miller, April 17, 2021
function secondary_destinations = Blend_waste(blend_progress_cells,SpecLimits,petro_props)
global PVPs VP_MWs
secondary_blending = [];
nspec = size(SpecLimits,2);
nfueltypes = size(SpecLimits,3);
blendstockVP_subset = PVPs(:,1:45);
blendstockMWs = VP_MWs(1:45);
blendstockGasVPs = PVPs(:,59:61);
blendstockDieselVPs = PVPs(:,62:63);
blendstockJetVPs = PVPs(:,64:end);

GasVPs = [blendstockVP_subset blendstockGasVPs];
DieselVPs = [blendstockVP_subset blendstockDieselVPs];
JetVPs = [blendstockVP_subset blendstockJetVPs];

blendstockGasMWs = VP_MWs(59:61);
blendstockDieselMWs = VP_MWs(62:63);
blendstockJetMWs = VP_MWs(64:end);

GasMWs = [blendstockMWs blendstockGasMWs];
DieselMWs = [blendstockMWs blendstockDieselMWs];
JetMWs = [blendstockMWs blendstockJetMWs];

for i=1:size(blend_progress_cells,1)
    secondary_mix = cell2mat(blend_progress_cells(i,2));
    if isempty(secondary_mix)
        secondary_blending = {"Diesel", "N/A", "N/A", "N/A", "N/A", "Jet", "N/A", "N/A", "N/A", "N/A", "Auto", "N/A", "N/A", "N/A", "N/A"};
    elseif sum(secondary_mix(:,3)) == 0 % If there's nothing in the mixture, don't mix it lol
        secondary_blending = {"Diesel", "N/A", "N/A", "N/A", "N/A", "Jet", "N/A", "N/A", "N/A", "N/A", "Auto", "N/A", "N/A", "N/A", "N/A"};
    else
        secondary_props0 = alcmix(secondary_mix,blendstockVP_subset);
        secondary_props_diesel = [secondary_props0(1:3) secondary_props0(9) secondary_props0(5:8)];
        secondary_props_jet = [secondary_props0(1:3) secondary_props0(10) secondary_props0(5:8)];
        secondary_props_gasoline = [secondary_props0(1:8)];
        secondary_meetspecs = zeros(nfueltypes,nspec);
        
        secondary_meetspecs(1,:) = checkspecs(secondary_props_diesel,SpecLimits(:,:,1));
        secondary_meetspecs(2,:) = checkspecs(secondary_props_jet,SpecLimits(:,:,2));
        secondary_meetspecs(3,:) = checkspecs(secondary_props_gasoline,SpecLimits(:,:,3));
%         for i=1:nfueltypes
%             secondary_meetspecs(i,:) = checkspecs(secondary_props,SpecLimits(:,:,i));
%         end
        
        if sum(secondary_meetspecs(1,:)) == nspec
            diesel_blend_max2 = 1;
            diesel_mix_props2 = secondary_props_diesel;
            dFailedProp2 = "Neat fuel";
            dHiLo2 = "Neat fuel";
        else
            [diesel_blend_max2, dmix_props2, dFailedProp2, dHiLo2] = Fuel_Blender(secondary_props_diesel, petro_props(1,:), SpecLimits(:,:,1), 0, "alcohol", "Diesel", secondary_mix);
            diesel_mix_props2 = dmix_props2(end,2:end);
        end
        
        if sum(secondary_meetspecs(2,:)) == nspec
            jet_blend_max2 = 1;
            jet_mix_props2 = secondary_props_jet;
            jFailedProp2 = "Neat fuel";
            jHiLo2 = "Neat fuel";
        else
            [jet_blend_max2, jmix_props2, jFailedProp2, jHiLo2] = Fuel_Blender(secondary_props_jet, petro_props(2,:), SpecLimits(:,:,2), 0, "alcohol", "Jet", secondary_mix);
            jet_mix_props2 = jmix_props2(end,2:end);
        end
        
        if sum(secondary_meetspecs(3,:)) == nspec
            auto_blend_max2 = 1;
            auto_mix_props2 = secondary_props_gasoline;
            aFailedProp2 = "Neat fuel";
            aHiLo2 = "Neat fuel";
        else
            [auto_blend_max2, amix_props2, aFailedProp2, aHiLo2] = Fuel_Blender(secondary_props_gasoline, petro_props(3,:), SpecLimits(:,:,3), 0, "alcohol", "Gasoline", secondary_mix);
            auto_mix_props2 = amix_props2(end,2:end);
        end
        
        secondary_blending_loop = {"Diesel", diesel_blend_max2, diesel_mix_props2, dFailedProp2, dHiLo2,...
            "Jet", jet_blend_max2, jet_mix_props2, jFailedProp2, jHiLo2,...
            "Auto", auto_blend_max2, auto_mix_props2, aFailedProp2, aHiLo2};
        secondary_blending = vertcat(secondary_blending, secondary_blending_loop);
    end
end
secondary_destinations = horzcat(blend_progress_cells, secondary_blending);
end