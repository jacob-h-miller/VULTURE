%  Make cuts from a biofuel mixture in order to mix it at higher weight
%  fractions with a petrofuel.
%  Jake Miller, April 17, 2021
function [blend_progress blend_progress_HCs] = Cut_Blend(VFAs, blend_max, FailedProp, HiLo, Cut_TopBot, HC_HiLo, HC_Cut_TopBot, HC_blend_max, HC_FailedProp, petro_props, Spec_Limits, fuel_class)

global PVPs

nspec = size(Spec_Limits,2);
working_alcs = VFAs;

alc_PVPs = PVPs(:,1:45);
HC_PVPs = PVPs(:,46:58);
if fuel_class == "Gasoline"
    fuel_PVPs = PVPs(:,59:61);
elseif fuel_class == "Diesel"
    fuel_PVPs = PVPs(:,62:63);
elseif fuel_class == "Jet"
    fuel_PVPs = PVPs(:,64:end);
else
    disp("Cut_Blend alcohol fuel class undefined");
end

alcPVPs = [alc_PVPs fuel_PVPs];
HCPVPs = [HC_PVPs fuel_PVPs];

alc_specs0 = alcmix(working_alcs, alc_PVPs);
if fuel_class == "Gasoline"
    alc_specs = alc_specs0(1:8);
elseif fuel_class == "Diesel"
    alc_specs = [alc_specs0(1:3) alc_specs0(9) alc_specs0(5:8)];
elseif fuel_class == "Jet"
    alc_specs = [alc_specs0(1:3) alc_specs0(10) alc_specs0(5:8)];
end
[~, working_specs0, ~, ~] = Fuel_Blender(alc_specs, petro_props, Spec_Limits, 0, "alcohol", fuel_class, working_alcs);
working_specs = working_specs0(end-1,2:end);
cut_alcs = working_alcs;
cut_alcs(:,3) = 0;
cut_alcs(:,13) = 0;
cut_alcs(:,end) = 0; % Make cut_alcs into a zero mass contents matrix that still has all the alcohol choices in it.
keep_cutting = "Yes";
blend_frac = blend_max;
blend_progress = {working_alcs, [], 1, blend_max, FailedProp, HiLo, working_specs, "N/A"};
if Cut_TopBot == "Bottom"
    alc_index = size(working_alcs,1);
elseif Cut_TopBot == "Top"
    alc_index = 1;
end
while keep_cutting == "Yes"
    % Cut either the top or bottom alcohol, and keep cutting if that
    % alcohol isn't present.
    if Cut_TopBot == "Bottom"
        cut_wtfrac = working_alcs(alc_index,end);
        while cut_wtfrac == 0
            %cut_alcs(alc_index,:) = working_alcs(alc_index,:);
            working_alcs(alc_index,3) = 0; %zero out the mole fraction,
            working_alcs(alc_index,13) = 0; % volume fraction,
            working_alcs(alc_index,end) = 0; % and mass fraction of the cut alcohol
            alc_index = alc_index-1;
            cut_wtfrac = working_alcs(alc_index,end); % Find weight fraction of next component to cut.
        end
        % This final cut will be of a component that's actually there.
        cut_alcs(alc_index,:) = working_alcs(alc_index,:);
        working_alcs(alc_index,3) = 0;
        working_alcs(alc_index,13) = 0;
        working_alcs(alc_index,end) = 0;
        alc_index = alc_index - 1;
    elseif Cut_TopBot == "Top"
        cut_wtfrac = working_alcs(alc_index,end);
        while cut_wtfrac == 0
            cut_alcs(alc_index,:) = working_alcs(alc_index,:);
            working_alcs(alc_index,3) = 0; %zero out the mole fraction,
            working_alcs(alc_index,13) = 0; % volume fraction,
            working_alcs(alc_index,end) = 0; % and mass fraction of the cut alcohol
            alc_index = alc_index+1;
            cut_wtfrac = working_alcs(alc_index,end); % Find weight fraction of next component to cut.
        end
        % This final cut will be of a component that's actually there.
        cut_alcs(alc_index,:) = working_alcs(alc_index,:);
        working_alcs(alc_index,3) = 0;
        working_alcs(alc_index,13) = 0;
        working_alcs(alc_index,end) = 0;
        alc_index = alc_index + 1;
    end
    % First, see if this new mixture is a neat fuel. If not, try to blend this new mixture.
    cut_specs0 = alcmix(working_alcs, alc_PVPs);
    if fuel_class == "Diesel"
        cut_specs = [cut_specs0(1:3) cut_specs0(9) cut_specs0(5:8)];
    elseif fuel_class == "Jet"
        cut_specs = [cut_specs0(1:3) cut_specs0(10) cut_specs0(5:8)];
    elseif fuel_class == "Gasoline"
        cut_specs = cut_specs0(1:8);
    end
    cut_meetspecs = zeros(size(Spec_Limits));
    cut_meetspecs = checkspecs(cut_specs,Spec_Limits);
    %     working_HCs = HYD_alc(working_alcs);
    %     HCcut_specs = alcmix(working_HCs);
    %     HCcut_meetspecs = zeros(size(Spec_Limits));
    %     HCcut_meetspecs = checkspecs(HCcut_specs,Spec_Limits);
    if sum(cut_meetspecs) == nspec
        working_wtfrac = sum(working_alcs(:,end)); % Mass fraction of original mixture still in blend
        progress_addition = {working_alcs, cut_alcs, working_wtfrac, 1, "Neat fuel", "N/A", cut_specs, "N/A"};
        blend_progress = vertcat(blend_progress, progress_addition);
        keep_cutting = "No";
    elseif sum(working_alcs(:,end)) == 0; % If there are no alcohols left, stop cutting.
        keep_cutting = "No";
    else
        working_specs0 = alcmix(working_alcs, alc_PVPs);
        if fuel_class == "Diesel"
            working_specs = [working_specs0(1:3) working_specs0(9) working_specs0(5:8)];
        elseif fuel_class == "Jet"
            working_specs = [working_specs0(1:3) working_specs0(10) working_specs0(5:8)];
        elseif fuel_class == "Gasoline"
            working_specs = working_specs0(1:8);
        end
        [blend_max_loop, mix_props_loop, FailedProp_loop, HiLo_loop] = Fuel_Blender(working_specs, petro_props, Spec_Limits, blend_frac, "alcohol", fuel_class, working_alcs);
        % Enter specs from this blending series into the progress cell
        % array. Take the last set of specs that *worked*
        mix_props_final = blendmix(blend_max_loop,working_specs,petro_props, "alcohol", fuel_class, working_alcs);
        mix_props_fail = blendmix(blend_max_loop+0.01,working_specs,petro_props, "alcohol", fuel_class, working_alcs);
        working_wtfrac = sum(working_alcs(:,end)); % Mass fraction of original mixture still in blend
        progress_addition = {working_alcs, cut_alcs, working_wtfrac, blend_max_loop, FailedProp_loop, HiLo_loop, mix_props_final, mix_props_fail};
        blend_progress = vertcat(blend_progress, progress_addition);
        % See if the next move to hit specs necessitates taking a cut from the
        % top or bottom of the mixture.
        Cut_TopBot_loop = Top_Bottom_Cut(FailedProp_loop, HiLo_loop);
        blend_frac = blend_max_loop;
        if Cut_TopBot_loop == Cut_TopBot % If we need to cut from the same end we were already, great, go for it.
            keep_cutting = "Yes";
        else % If we now need to cut from a different end than we did already, stop. Only one cut total is allowed and we would now need two.
            keep_cutting = "No";
        end
        if blend_max_loop == 1 % If we end up with a fuel that's 100% blendable, return it!
            keep_cutting = "No";
        end
    end
end

blend_progress_HCs = {};
HC_blend_frac = 0;

if HC_Cut_TopBot == "Bottom"
    HC_alc_index = size(working_alcs,1);
elseif HC_Cut_TopBot == "Top"
    HC_alc_index = 1;
end

% for i=1:size(blend_progress,1)
%     working_HCs = HYD_alc(blend_progress{i,1});
%     working_specs_HCs = alcmix(working_HCs);
%     cut_alcs_HCs = blend_progress{i,2};
%     working_wtfrac_HCs = blend_progress{i,3};
%     [blend_max_loop_HCs, mix_props_loop_HCs, FailedProp_loop_HCs, HiLo_loop_HCs] = Fuel_Blender(working_specs_HCs, petro_props, Spec_Limits, HC_blend_frac);
%     mix_props_final_HCs = blendmix(blend_max_loop_HCs,working_specs_HCs,petro_props);
%     mix_props_fail_HCs = blendmix(blend_max_loop_HCs+0.01,working_specs_HCs,petro_props);
%     progress_addition_HCs = {working_HCs, cut_alcs_HCs, working_wtfrac_HCs, blend_max_loop_HCs, FailedProp_loop_HCs, HiLo_loop_HCs, mix_props_final_HCs, mix_props_fail_HCs};
%     blend_progress_HCs = vertcat(blend_progress_HCs, progress_addition_HCs);
% end

% Now, keep blending hydrogenated streams until we hit a criterion to stop.
HC_keep_cutting = "Yes";
working_alcs_HCs = VFAs;
cut_alcs = VFAs;
cut_alcs(:,3) = 0;
cut_alcs(:,13) = 0;
cut_alcs(:,end) = 0; % Make cut_alcs into a zero mass contents matrix that still has all the alcohol choices in it.

working_HCs = HYD_alc(working_alcs_HCs);
%blend_progress = {working_alcs, [], 1, blend_max, FailedProp, HiLo, working_specs, "N/A"};
alk_specs0 = HCmix(working_HCs,HC_PVPs);

if fuel_class == "Gasoline"
    alk_specs = alk_specs0(1:8);
elseif fuel_class == "Diesel"
    alk_specs = [alk_specs0(1:3) alk_specs0(9) alk_specs0(5:8)];
elseif fuel_class == "Jet"
    alk_specs = [alk_specs0(1:3) alk_specs0(10) alk_specs0(5:8)];
end

[~, working_specs0, ~, ~] = Fuel_Blender(alk_specs, petro_props, Spec_Limits, 0, "hydrocarbon", fuel_class, working_HCs);
working_specs = working_specs0(end-1,2:end);
blend_progress_HCs = {working_HCs, cut_alcs, 1, HC_blend_max, HC_FailedProp, HC_HiLo, working_specs, "N/A"};
while HC_keep_cutting == "Yes"
    % Cut either the top or bottom alcohol, and keep cutting if that
    % alcohol isn't present.
    if HC_Cut_TopBot == "Bottom" % still use the same cut location as for alcohols
        cut_wtfrac = working_alcs_HCs(HC_alc_index,end);
        while cut_wtfrac == 0
            %cut_alcs(HC_alc_index,:) = working_alcs_HCs(HC_alc_index,:);
            working_alcs_HCs(HC_alc_index,3) = 0;
            working_alcs_HCs(HC_alc_index,13) = 0;
            working_alcs_HCs(HC_alc_index,end) = 0;
            HC_alc_index = HC_alc_index - 1;
            cut_wtfrac = working_alcs_HCs(HC_alc_index,end); % Find weight fraction of next component to cut.
        end
        % This final cut will be of a component that's actually there.
        cut_alcs(HC_alc_index,:) = working_alcs_HCs(HC_alc_index,:);
        working_alcs_HCs(HC_alc_index,3) = 0;
        working_alcs_HCs(HC_alc_index,13) = 0;
        working_alcs_HCs(HC_alc_index,end) = 0;
        HC_alc_index = HC_alc_index - 1;
    elseif HC_Cut_TopBot == "Top"
        cut_wtfrac = working_alcs_HCs(HC_alc_index,end);
        while cut_wtfrac == 0
            if working_alcs_HCs(1,end) > 0
                disp('why tho');
            end
            cut_alcs(HC_alc_index,:) = working_alcs_HCs(HC_alc_index,:);
            working_alcs_HCs(HC_alc_index,3) = 0;
            working_alcs_HCs(HC_alc_index,13) = 0;
            working_alcs_HCs(HC_alc_index,end) = 0;
            HC_alc_index = HC_alc_index + 1;
            if HC_alc_index > size(working_alcs_HCs,1)
                disp('what the fuck');
            end
            cut_wtfrac = working_alcs_HCs(HC_alc_index,end); % Find weight fraction of next component to cut.
        end
        % This final cut will be of a component that's actually there.
        cut_alcs(HC_alc_index,:) = working_alcs_HCs(HC_alc_index,:);
        working_alcs_HCs(HC_alc_index,3) = 0;
        working_alcs_HCs(HC_alc_index,13) = 0;
        working_alcs_HCs(HC_alc_index,end) = 0;
        HC_alc_index = HC_alc_index + 1;
    end
    % First, see if this new mixture is a neat fuel. If not, try to blend this new mixture.
    if sum(working_alcs_HCs(:,end)) == 0 % If there isn't anything left to cut, stop cutting and leave the function.
        HC_keep_cutting = "No";
        return
    end
    working_HCs = HYD_alc(working_alcs_HCs);
    cut_specs0 = HCmix(working_HCs,HC_PVPs);
    if fuel_class == "Diesel"
        cut_specs = [cut_specs0(1:3) cut_specs0(9) cut_specs0(5:8)];
    elseif fuel_class == "Jet"
        cut_specs = [cut_specs0(1:3) cut_specs0(10) cut_specs0(5:8)];
    elseif fuel_class == "Gasoline"
        cut_specs = cut_specs0(1:8);
    end
    cut_meetspecs = zeros(size(Spec_Limits));
    cut_meetspecs = checkspecs(cut_specs,Spec_Limits);
    blend_frac = 0;
    if sum(cut_meetspecs) == nspec
        working_wtfrac = sum(working_alcs_HCs(:,end)); % Mass fraction of original mixture still in blend
        progress_addition_HCs = {working_HCs, cut_alcs, working_wtfrac, 1, "Neat fuel", "N/A", cut_specs, "N/A"};
        blend_progress_HCs = vertcat(blend_progress_HCs, progress_addition_HCs);
        HC_keep_cutting = "No";
    else
        working_specs0 = alcmix(working_HCs, HC_PVPs);
        if fuel_class == "Diesel"
            working_specs = [working_specs0(1:3) working_specs0(9) working_specs0(5:8)];
        elseif fuel_class == "Jet"
            working_specs = [working_specs0(1:3) working_specs0(10) working_specs0(5:8)];
        elseif fuel_class == "Gasoline"
            working_specs = working_specs0(1:8);
        end
        [blend_max_loop, mix_props_loop, FailedProp_loop, HiLo_loop] = Fuel_Blender(working_specs, petro_props, Spec_Limits, blend_frac, "hydrocarbon", fuel_class, working_HCs);
        % Enter specs from this blending series into the progress cell
        % array. Take the last set of specs that *worked*
        mix_props_final = blendmix(blend_max_loop,working_specs,petro_props, 'hydrocarbon', fuel_class, working_HCs);
        mix_props_fail = blendmix(blend_max_loop+0.01,working_specs,petro_props, 'hydrocarbon', fuel_class, working_HCs);
        working_wtfrac = sum(working_alcs_HCs(:,end)); % Mass fraction of original mixture still in blend
        if blend_max_loop > 1
            disp(">100% fuel");
        end
        progress_addition = {working_HCs, cut_alcs, working_wtfrac, blend_max_loop, FailedProp_loop, HiLo_loop, mix_props_final, mix_props_fail};
        blend_progress_HCs = vertcat(blend_progress_HCs, progress_addition);
        % See if the next move to hit specs necessitates taking a cut from the
        % top or bottom of the mixture.
        Cut_TopBot_loop = Top_Bottom_Cut(FailedProp_loop, HiLo_loop);
        blend_frac = blend_max_loop;
        if Cut_TopBot_loop == HC_Cut_TopBot % If we need to cut from the same end we were already, great, go for it.
            HC_keep_cutting = "Yes";
        else % If we now need to cut from a different end than we did already, stop. Only one cut total is allowed and we would now need two.
            HC_keep_cutting = "No";
        end
        if blend_max_loop == 1 % If we end up with a fuel that's 100% blendable, return it!
            HC_keep_cutting = "No";
        end
    end
end