% Logic for fuel uses for various volatile fatty acid-derived alcohol
% mixtures, with separation following reaction. Jake Miller, 5/3/2021

function [VFAmatrix, inlet_specs_alc0, dmix_props, jmix_props, amix_props, diesel_blend_max, jet_blend_max, auto_blend_max] = FuelCriteria_noseps(VFA_molfracs) 

global PVPs VP_MWs

%tic
%Create a vector of VFA alcohol concentration
%VFA species in order in VFA_molfracs: acetic, propionic, isobutyric,
%butyric, isopentanoic, pentanoic, hexanoic, heptanoic, ocatanoic,
%nonanoic, decanoic
%VFA_molfracs = [1 0 0 0 0 0 0 0 0 0 0];

[VFACarbon, VFA_old, Light,Heavy]= KetoneGenDriver(VFA_molfracs);

% Take the list of alcohols from KetoneGenDriver (45 alcohols, including
% one duplicate, no discernable order) and turn it into a new list (44
% alcohols, no duplicates, ordered by boiling point).
VFA = VFAreorder(VFA_old)'; 

VFAComp = VFA/sum(VFA);
nalcohols = size(VFA,1);
%%Read data from the idividual predicted fuel property spreadsheet
% database = readtable('Prediction_FP_Database.xlsx', 'Range', 'A1:BT48');
% opts = detectImportOptions('Prediction_FP_Database.xlsx');
% preview('Prediction_FP_Database.xlsx',opts);

%List of VFA alcohols

VFAalcohols= readcell('Prediction_FP_Database.xlsx', 'Range', 'A4:A48');
VFAindex=linspace(1,length(VFAalcohols),length(VFAalcohols))';
%read molecular weights /(g/mol)
MWs = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'Q4:Q48');
%read flashpoint (degree C)
FP = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'AX4:AX48');
%read energy density (LHV), MJ/kg
LHV = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'AN4:AN48');
%read energy density (HHV), MJ/kg
HHV = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'Af4:Af48');
%read boiling point, C
BP = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'r4:r48');
%read viscosity, cP
U = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'bb4:bb48');
% read melting point, C
MP = readmatrix('Prediction_FP_Database.xlsx', 'Range', 't4:t48');
% read vapor pressure at 25C, in mmhg
VP = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'v4:v48');
% read water solubility, in mg/L
Wsol = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'x4:x48');
% read yield sooting index
YSI = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'z4:z48');
% read mass% H
Hpct = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'ai4:ai48');
% read density /(g/mL)
rho = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'az4:az48');
% read cetane number
CN = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'bd4:bd48');
% read research octane number
RON = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'bf4:bf48');
% read motor octane number
MON = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'bi4:bi48');
%read the carbon number of each VFA alcohol
CarbonNumber = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'G4:G48');
%read the Torabian and Sobati blending index
BI_general = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'BU4:BU48');
%read viscosities at 40 and -40 C
U40 = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'BV4:BV48');
Uneg40 = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'BW4:BW48');
wtfracs = moltowtfrac(VFAComp,MWs);
volfracs = volfraccalc(wtfracs,rho);

% mass fractions of linear alcohols
wtfracs_forplot_linear = zeros(13,2);
wtfracs_forplot_linear(:,1) = 3:15;
wtfracs_forplot_linear(1,2) = wtfracs(1);
wtfracs_forplot_linear(2,2) = wtfracs(2);
wtfracs_forplot_linear(3,2) = wtfracs(4) + wtfracs(5);
wtfracs_forplot_linear(4,2) = wtfracs(9) + wtfracs(10);
wtfracs_forplot_linear(5,2) = sum(wtfracs(14:16));
wtfracs_forplot_linear(6,2) = sum(wtfracs(20:22));
wtfracs_forplot_linear(7,2) = sum(wtfracs(25:28));
wtfracs_forplot_linear(8,2) = sum(wtfracs(31:33));
wtfracs_forplot_linear(9,2) = sum(wtfracs(36:38));
wtfracs_forplot_linear(10,2)= sum(wtfracs(40:41));
wtfracs_forplot_linear(11,2)= sum(wtfracs(42:43));
wtfracs_forplot_linear(12,2)= wtfracs(44);
wtfracs_forplot_linear(13,2)= wtfracs(45);

%mass fractions of branched alcohols
wtfracs_forplot_branched = zeros(13,2);

wtfracs_forplot_branched(:,1) = 3:15;
wtfracs_forplot_branched(3,2) = wtfracs(3);
wtfracs_forplot_branched(4,2) = sum(wtfracs(6:7));
wtfracs_forplot_branched(5,2) = wtfracs(8) + sum(wtfracs(11:12));
wtfracs_forplot_branched(6,2) = wtfracs(13) + sum(wtfracs(17:18));
wtfracs_forplot_branched(7,2) = wtfracs(19) + sum(wtfracs(23:24));
wtfracs_forplot_branched(8,2) = sum(wtfracs(29:30));
wtfracs_forplot_branched(9,2) = sum(wtfracs(34:35));
wtfracs_forplot_branched(10,2)= wtfracs(39);

% figure
% hold on
% bar(wtfracs_forplot_linear(:,1), 100*[wtfracs_forplot_linear(:,2) wtfracs_forplot_branched(:,2)],'stacked');
% cnumber_wtpct = (wtfracs_forplot_linear(:,2) + wtfracs_forplot_branched(:,2))*100;
% xlim([2.5 15.5]);
% Cnumbers = num2cell(wtfracs_forplot_linear(:,1));
% labels = arrayfun(@(value) num2str(value,2),cnumber_wtpct,'UniformOutput',false);
% text(wtfracs_forplot_linear(:,1),cnumber_wtpct,labels,...
%   'HorizontalAlignment','center',...
%   'VerticalAlignment','bottom');
% xlabel('Alcohol Carbon Number');
% ylabel('Mass Percent');
% legend('Linear alcohols','Branched alcohols','Location','Northwest');

%VFAmatrix_full =[VFAindex, CarbonNumber, VFAComp, BP, FP, LHV, U, MP, VP, Wsol, YSI, rho, CN, RON, wtfracs];
VFAmatrix = [VFAindex, CarbonNumber, VFAComp, BP, FP, LHV, U, MP, Wsol, CN, RON, rho, volfracs, BI_general, U40, Uneg40, wtfracs];
% Lower (first row) and upper (second row) limits for flash point (deg C),
% LHV (MJ/kg), boiling point (deg C), viscosity (cP).
SpecLimits = zeros(2,8,3); % Dimensions: 2-lower and upper; 8-number of properties; 3-number of fuel types
% specs need to be in the same order as columns in VFA matrix. If there are
% no upper/lower limits, putting -1000 and 1000, except for water
% solubility, the upper limit of which is 1000000 mg/L (that is, 100%).
SpecLimits(1,:,1) = [-1000, 52, 36, 1.615, -1000, -1000, 40, -1000]; % Diesel lower limits
SpecLimits(2,:,1) = [338, 1000, 1000, 3.485, 0, 20000, 1000, 1000]; % Diesel upper limits
SpecLimits(1,:,2) = [-1000, 38, 42, -1000, -1000, -1000, 45, -1000]; % Jet lower limits
SpecLimits(2,:,2) = [300, 1000, 1000, 10.08, -40, 1000000, 1000, 1000]; % Jet upper limits
SpecLimits(1,:,3) = [-1000, -1000, 26.8, -1000, -1000, -1000, -1000, 87.5]; % Auto lower limits
SpecLimits(2,:,3) = [171, 1000, 1000, 1000, -10, 1000000, 1000, 1000]; % Auto upper limits
nspec = size(SpecLimits,2);
nfueltypes = size(SpecLimits,3);
comp_specs = zeros(nalcohols,nspec,nfueltypes); % Each entry indicates whether or not the alcohol meets each spec by itself.
% Enumerate which fuel components meet which specs
for i=1:nalcohols
    for j=1:nspec
        VFAprop = VFAmatrix(i,3+j);
        for k=1:nfueltypes
            if (VFAprop > SpecLimits(1,j,k)) && (VFAprop < SpecLimits(2,j,k))
                comp_specs(i,j,k) = 1;
            end
        end
    end
end

% Properties of initial fuel mixture. Calculate FP and viscosity with
% carbon number blend properties, all else linearly, using the alcmix
% function.

alc_PVPs = PVPs(:,1:45);
diesel_PVPs = PVPs(:,62:63);
gasoline_PVPs = PVPs(:,59:61);
jet_PVPs = PVPs(:,64:end);

dmix_PVPs = [alc_PVPs diesel_PVPs];
gmix_PVPs = [alc_PVPs gasoline_PVPs];
jmix_PVPs = [alc_PVPs jet_PVPs];

inlet_specs_alc0 = alcmix(VFAmatrix,alc_PVPs);
% customize inlet specs with appropriate viscosity (40C for diesel, -40C
% for jet, and room temp (not super relevant) for gasoline
inlet_specs_alc_diesel = [inlet_specs_alc0(1:3) inlet_specs_alc0(9) inlet_specs_alc0(5:8)];
inlet_specs_alc_jet = [inlet_specs_alc0(1:3) inlet_specs_alc0(10) inlet_specs_alc0(5:8)];
inlet_specs_alc_gasoline = [inlet_specs_alc0(1:8)];

% 1) Determine which specs the inlet fuel meets for each of the three fuel
% types. If any of them meet all the specs, indicate that through
% Neat_#fueltype = 1.

Neat_diesel = 0;
Neat_jet = 0;
Neat_auto = 0;

secondary_meetspecs = zeros(nfueltypes,nspec);

secondary_meetspecs(1,:) = checkspecs(inlet_specs_alc_diesel,SpecLimits(:,:,1));
secondary_meetspecs(2,:) = checkspecs(inlet_specs_alc_diesel,SpecLimits(:,:,2));
secondary_meetspecs(3,:) = checkspecs(inlet_specs_alc_diesel,SpecLimits(:,:,3));

if sum(secondary_meetspecs(1,:)) == nspec
    Neat_diesel = 1;
end

if sum(secondary_meetspecs(2,:)) == nspec
    Neat_jet = 1;
end

if sum(secondary_meetspecs(3,:)) == nspec
    Neat_auto = 1;
end

% in petro properties, properties are listed in the order: diesel, jet,
% gasoline

petro_BP = readmatrix('Neat_Petrofuel_Properties.xlsx', 'Range', 'C2:C4');
petro_FP = readmatrix('Neat_Petrofuel_Properties.xlsx', 'Range', 'AG2:AG4');
petro_LHV = readmatrix('Neat_Petrofuel_Properties.xlsx', 'Range', 'Y2:Y4');
petro_U = readmatrix('Neat_Petrofuel_Properties.xlsx', 'Range', 'AM2:AM4');
petro_MP = readmatrix('Neat_Petrofuel_Properties.xlsx', 'Range', 'E2:E4');
petro_Wsol = readmatrix('Neat_Petrofuel_Properties.xlsx', 'Range', 'I2:I4');
petro_CN = readmatrix('Neat_Petrofuel_Properties.xlsx', 'Range', 'AO2:AO4');
petro_RON = readmatrix('Neat_Petrofuel_Properties.xlsx','Range', 'AQ2:AQ4');

petro_props = [petro_BP, petro_FP, petro_LHV, petro_U, petro_MP, petro_Wsol, petro_CN, petro_RON];
%prop_list = ["Boiling Point","Flash Point","Lower Heating Value","Viscosity","Melting Point","Water Solubility","Cetane Number"];

% 2) Start blending! See if fuel can meet specs with 10% (mass), then keep
% adding on 1%s until it doesn't meet specs. Assume linear blending by
% weight percent for all specs, since we don't have carbon numbers for
% example petrofuels.

% First, try to blend into diesel.
if Neat_diesel == 1
    diesel_blend_max = 1;
    dmix_props = 'Neat fuel';
else
    start_frac = 0;
    [diesel_blend_max, dmix_props, dFailedProp, dHiLo] = Fuel_Blender(inlet_specs_alc_diesel, petro_props(1,:), SpecLimits(:,:,1), start_frac, "alcohol", "Diesel", VFAmatrix);
    Diesel_cut = Top_Bottom_Cut(dFailedProp, dHiLo);
end

% Then, try to blend into jet fuel.
if Neat_jet == 1
    jet_blend_max = 1;
    jmix_props = 'Neat fuel';
else
    start_frac = 0;
    [jet_blend_max, jmix_props, jFailedProp, jHiLo] = Fuel_Blender(inlet_specs_alc_jet, petro_props(2,:), SpecLimits(:,:,2), start_frac, "alcohol", "Jet", VFAmatrix);
    Jet_cut = Top_Bottom_Cut(jFailedProp, jHiLo);
end

% Finally, try to blend into auto fuel.
if Neat_auto == 1
    auto_blend_max = 1;
    amix_props = 'Neat fuel';
else
    start_frac = 0;
    [auto_blend_max, amix_props, aFailedProp, aHiLo] = Fuel_Blender(inlet_specs_alc_gasoline, petro_props(3,:), SpecLimits(:,:,3), start_frac, "alcohol", "Gasoline", VFAmatrix);
    Auto_cut = Top_Bottom_Cut(aFailedProp, aHiLo);
end

%toc
end
 
