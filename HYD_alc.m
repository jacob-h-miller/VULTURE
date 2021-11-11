% Takes alcohols and hydrogenates them to alkanes; associates each species
% with relevant fuel and physiochemical properties.

function HYD_results = HYD_alc(alcohols);

VFAalkanes= readcell('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'b4:b16');
Alkaneindex=linspace(1,length(VFAalkanes),length(VFAalkanes))';
%read molecular weights /(g/mol)
MWs = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'r4:r16');
%read flashpoint (degree C). From UMass Lowell spreadsheet or EPA TEST
FP = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'bf4:bf16');
%read energy density (LHV), MJ/kg. Average of LD (%H), Dulong, and Boie
%methods.
LHV = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'Aq4:Aq16');
%read energy density (HHV), MJ/kg
HHV = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'Ai4:Ai16');
%read boiling point, C. From EPISuite/EPA TEST
BP = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 't4:t16');
%read viscosity, cP. From EPA TEST
U = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'bo4:bo16');
% read melting point, C. From EPISuite/EPA TEST
MP = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'v4:v16');
% read vapor pressure at 25C, in mmhg
VP = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'x4:x16');
% read water solubility, in mg/L. From EPISuite/EPA TEST
Wsol = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'z4:z16');
% read yield sooting index
YSI = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'aa4:aa16');
% read mass% H
Hpct = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'al4:al16');
% read density /(g/mL)
rho = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'bg4:bg16');
% read cetane number. From UMass Lowell database.
CN = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'bt4:bt16');
% read research octane number. From UMass Lowell database.
RON = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'bv4:bv16');
% read motor octane number. From UMass Lowell database.
MON = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'bx4:bx16');
%read the carbon number of each VFA alcohol.
CarbonNumber = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'h4:h16');
%read the Wickey-Chittenden blending indices
BI_HCs = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'CG4:CG16');
%read the Torabian and Sobati blending indices
BI_general = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'CI4:CI16');
%read the viscosities at -40C and 40C.
U40 = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'CK4:CK16');
Uneg40 = readmatrix('Biorefinery_hydrocarbon_database.xlsx', 'Range', 'CJ4:CJ16');
molfracs = zeros(length(MWs),1);

for i=1:size(alcohols,1)
    Cnumber = alcohols(i,2);
    mol_frac = alcohols(i,3);
    molfracs(Cnumber-2) = molfracs(Cnumber-2) + mol_frac;
end

wtfracs = moltowtfrac(molfracs,MWs);
volfracs = volfraccalc(wtfracs, rho);
% Note: In the current program, we're never going to blend these
% hydrocarbons with oxygenates, so don't need to carry the general blending
% index number in the component vector.
HYD_results = [Alkaneindex, CarbonNumber, molfracs, BP, FP, LHV, U, MP, Wsol, CN, RON, rho, volfracs, BI_HCs, U40, Uneg40, wtfracs];

end
