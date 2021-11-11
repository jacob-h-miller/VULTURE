% Determine specs of a blended mix, where we specify properties of the
% alcohol mixture, the petrofuel, and the blend weight fractions.
%  Jake Miller, April 13, 2021. Last edited June 21, 2021.
function blendmixspecs = blendmix(blend_frac, alcprops, petroprops, blendstock_type, fuel_class, biomatrix)
global PVPs VP_MWs

if size(alcprops,2) > 8
    disp('alcprops too big');
end

nprops = size(alcprops,2);
blend_mat = zeros(2,nprops+1); % first row of blend_mat is props of alcohol fuel, second row is props of petrofuel.
blend_mat(:,1) = [blend_frac; 1-blend_frac];
blend_mat(:,2:end) = [alcprops; petroprops];
blendmixspecs = zeros(1,nprops);
if blendstock_type == "alcohol"
    blendstockVP_subset = PVPs(:,1:45);
    blendstockMWs = VP_MWs(1:45);
elseif blendstock_type == "hydrocarbon"
    blendstockVP_subset = PVPs(:,46:58);
    blendstockMWs = VP_MWs(46:58);
else
    disp('Blendstock type for vapor pressure calculation is unspecified.');
end

if fuel_class == "Gasoline"
    fuelVP_subset = PVPs(:,59:61);
    fuelMWs = VP_MWs(59:61);
    fuelwtfracs = [0.25; 0.5; 0.25]; % toluene, trimethylpentane, and propylbenzene
    fueldensity = [0.744; 0.744; 0.744]; % To make life easier, just use total fuel density for everything
    fuelU = [petroprops(4); petroprops(4); petroprops(4)];
elseif fuel_class == "Diesel"
    fuelVP_subset = PVPs(:,62:63);
    fuelMWs = VP_MWs(62:63);
    fuelwtfracs = [0.68; 0.32]; % nonadecane and eicosane
    fueldensity = [0.853; 0.853]; % To make life easier, just use total fuel density for everything
    fuelU = [petroprops(4); petroprops(4)];
elseif fuel_class == "Jet"
    fuelVP_subset = PVPs(:,64:end);
    fuelMWs = VP_MWs(64:end);
    fuelwtfracs = [1]; % Single "jet fuel" component
    fueldensity = [0.802];
    fuelU = [petroprops(4)];
else
    disp('Petrofuel type for vapor pressure calculation is unspecified');
end
% Matrix of vapor pressures as a function of temperature to be used in RaoultsLawCalc for boiling point
VP_matrix = [blendstockVP_subset fuelVP_subset];
% Vector of molecular weights as a function of all of these compounds
VP_molwts = [blendstockMWs fuelMWs];
biowtfracs = biomatrix(:,end);
blend_biowtfracs = blend_frac * biowtfracs;
blend_fuelwtfracs = (1-blend_frac) * fuelwtfracs;
blend_wtfracs = [blend_biowtfracs; blend_fuelwtfracs]; % Weight fraction of all "molecules" in blended fuel
blend_molfracs = wttomolfrac(blend_wtfracs, VP_molwts'); % Mole fraction of all "molecules" in blended fuel
alcohol_BIs = [2.5626E-103 1.2869E-102 2.8773E-105]; % Blending indices for diesel, jet A, and gasoline (respectively) with alcohols (formula from Torabian and Sobati)
HC_BIs = [0.014358077 0.018108274 0.007599803]; % Blending indices for diesel, jet A, and gasoline (respectively) with hydrocarbons (formula from Wickey and Chittenden)
if blendstock_type == "alcohol"
    if fuel_class == "Gasoline"
        fuel_BI = alcohol_BIs(3);
    elseif fuel_class == "Diesel"
        fuel_BI = alcohol_BIs(1);
    elseif fuel_class == "Jet"
        fuel_BI = alcohol_BIs(2);
    else
        disp("Blend index: Petrofuel class not specified.");
    end
elseif blendstock_type == "hydrocarbon"
    if fuel_class == "Gasoline"
        fuel_BI = HC_BIs(3);
    elseif fuel_class == "Diesel"
        fuel_BI = HC_BIs(1);
    elseif fuel_class == "Jet"
        fuel_BI = HC_BIs(2);
    else
        disp("Blend index: Petrofuel class not specified.");
    end
else
    disp("Blend index: Blendstock type not specified.");
end
    
% First, calculate boiling point
blendmixspecs(1) = RaoultsLawCalc(blend_molfracs,VP_matrix);
% Then, calculate flash point. Use the Torabian and Sobati Blending Index formula if the
% blendstock is an alcohol and use the Wickey and Chittenden Blending Index
% formula if the blendstock is a hydrocarbon.
blendBIs = zeros(size(blend_molfracs));
blendBIs(1:size(biomatrix,1)) = biomatrix(:,14);
blendBIs(size(biomatrix,1)+1:end) = fuel_BI;

if blendstock_type == "alcohol"
    blendmixspecs(2) = TorabianSobatiFP(blend_molfracs,blendBIs);
elseif blendstock_type == "hydrocarbon"
    blendrhos = zeros(size(blend_molfracs));
    blendrhos = [biomatrix(:,12); fueldensity];
    blendvolfracs = volfraccalc(blend_wtfracs, blendrhos);
    blendmixspecs(2) = WickeyChittendenFP(blendvolfracs, blendBIs);
else
    disp("Blendstock type for flash point tabulation not specified");
end

% Then, calculate lower heating value
blendmixspecs(3) = LinearBlend(blend_mat(:,4), blend_mat(:,1));

% Then, calculate viscosity
if strcmp(fuel_class,'Gasoline')
    bioU = biomatrix(:,7); % viscosity at ambient conditions
elseif strcmp(fuel_class,'Diesel')
    bioU = biomatrix(:,15); % viscosity at 40C
elseif strcmp(fuel_class,'Jet')
    bioU = biomatrix(:,16); % viscosity at -40C
end
blendUs = [bioU; fuelU];

blendmixspecs(4) = GNviscosity(blend_molfracs,blendUs);

% Then, calculate melting point, water solubility, CN, and RON (all linear
% blends)

for i=5:8
    blendmixspecs(i) = LinearBlend(blend_mat(:,i+1), blend_mat(:,1));
end

end


