% Main function for VFA upgrading scheme. Can either react, then separate
% (FuelCriteria_postrxnseps.m) or can separate, then react. Jake Miller,
% 5/3/2021.

function [VFAmatrix, inlet_specs_alc, inlet_specs_HC,  diesel_blend_secondary, jet_blend_secondary, auto_blend_secondary, dmix_props, jmix_props, amix_props, noseps_props] = FuelCriteriaMain()
clear all
close all
tic

global PVPs VP_MWs

% Order: C2, C3, i-C4, C4, i-C5, C5, C6, C7, C8
VFA_molfracs = [0.0082 0.016 0 0.239 0 0.012 0.725 0 0];
VFA_MWs = [60.05 74.08 88.11 88.11 102.13 102.13 116.16 130.18 144.21]; 

PVPs = readmatrix('Molecule_vapor_pressures.xlsx', 'Range', 'C10:BO452');
VP_MWs = readmatrix('VPs_MWs.xlsx', 'Range', 'A2:BM2');

[VFAmatrix, inlet_specs_alc, inlet_specs_HC, diesel_blend_secondary, jet_blend_secondary, auto_blend_secondary, dmix_props, jmix_props, amix_props] = FuelCriteria_postrxnseps(VFA_molfracs);

%read molecular weights from alcohol database /(g/mol)
MWs = readmatrix('Prediction_FP_Database.xlsx', 'Range', 'Q4:Q48');
nVFAs = size(VFA_molfracs,2);
noseps_props = cell(nVFAs-1,34);

for i=1:nVFAs-1
    VFA_molfracs_small = zeros(size(VFA_molfracs));
    VFA_molfracs_large = zeros(size(VFA_molfracs));
    VFA_molfracs_small(1:i) = VFA_molfracs(1:i);
    VFA_molfracs_large(1+i:end) = VFA_molfracs(1+i:end);
    small_mol = sum(VFA_molfracs_small)/2;
    large_mol = sum(VFA_molfracs_large)/2;
    if (small_mol == 0)
        [VFAmatrixsmall, inlet_specssmall, dmix_propssmall, jmix_propssmall, amix_propssmall, diesel_blend_maxsmall, jet_blend_maxsmall, auto_blend_maxsmall] = deal('N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A');
        [VFAmatrixsmall_HC, inlet_specssmall_HC, dmix_propssmall_HC, jmix_propssmall_HC, amix_propssmall_HC, diesel_blend_maxsmall_HC, jet_blend_maxsmall_HC, auto_blend_maxsmall_HC] = deal('N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A');
    else
        [VFAmatrixsmall, inlet_specssmall, dmix_propssmall, jmix_propssmall, amix_propssmall, diesel_blend_maxsmall, jet_blend_maxsmall, auto_blend_maxsmall] = FuelCriteria_noseps(VFA_molfracs_small);
        [VFAmatrixsmall_HC, inlet_specssmall_HC, dmix_propssmall_HC, jmix_propssmall_HC, amix_propssmall_HC, diesel_blend_maxsmall_HC, jet_blend_maxsmall_HC, auto_blend_maxsmall_HC] = FuelCriteria_noseps_HYD(VFA_molfracs_small);
    end
    noseps_props_small = {VFAmatrixsmall, inlet_specssmall, dmix_propssmall, jmix_propssmall, amix_propssmall, diesel_blend_maxsmall, jet_blend_maxsmall, auto_blend_maxsmall};
    noseps_props_small_HC = {VFAmatrixsmall_HC, inlet_specssmall_HC, dmix_propssmall_HC, jmix_propssmall_HC, amix_propssmall_HC, diesel_blend_maxsmall_HC, jet_blend_maxsmall_HC, auto_blend_maxsmall_HC};
    if (large_mol == 0)
        [VFAmatrixlarge, inlet_specslarge, dmix_propslarge, jmix_propslarge, amix_propslarge, diesel_blend_maxlarge, jet_blend_maxlarge, auto_blend_maxlarge] = deal('N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A');
        [VFAmatrixlarge_HC, inlet_specslarge_HC, dmix_propslarge_HC, jmix_propslarge_HC, amix_propslarge_HC, diesel_blend_maxlarge_HC, jet_blend_maxlarge_HC, auto_blend_maxlarge_HC] = deal('N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A');
    else
        [VFAmatrixlarge, inlet_specslarge, dmix_propslarge, jmix_propslarge, amix_propslarge, diesel_blend_maxlarge, jet_blend_maxlarge, auto_blend_maxlarge] = FuelCriteria_noseps(VFA_molfracs_large);
        [VFAmatrixlarge_HC, inlet_specslarge_HC, dmix_propslarge_HC, jmix_propslarge_HC, amix_propslarge_HC, diesel_blend_maxlarge_HC, jet_blend_maxlarge_HC, auto_blend_maxlarge_HC] = FuelCriteria_noseps_HYD(VFA_molfracs_large);
    end
    noseps_props_large = {VFAmatrixlarge, inlet_specslarge, dmix_propslarge, jmix_propslarge, amix_propslarge, diesel_blend_maxlarge, jet_blend_maxlarge, auto_blend_maxlarge};
    noseps_props_large_HC = {VFAmatrixlarge_HC, inlet_specslarge_HC, dmix_propslarge_HC, jmix_propslarge_HC, amix_propslarge_HC, diesel_blend_maxlarge_HC, jet_blend_maxlarge_HC, auto_blend_maxlarge_HC};
    noseps_props(i,1:8) = noseps_props_small;
    noseps_props(i,9:16) = noseps_props_small_HC;
    noseps_props(i,18:25) = noseps_props_large;
    noseps_props(i,26:33) = noseps_props_large_HC;
    if small_mol == 0
        small_mass = 0;
    else
        small_mass = small_mol * dot(MWs,VFAmatrixsmall(:,3));
    end
    if large_mol == 0
        large_mass = 0;
    else
        large_mass = large_mol * dot(MWs,VFAmatrixlarge(:,3));
    end
    if isnan(small_mass)
        small_mass = 0;
    end
    
    if isnan(large_mass)
        large_mass = 0;
    end
    smallwtfrac = small_mass / (small_mass + large_mass);
    largewtfrac = large_mass / (small_mass + large_mass);
    
    noseps_props(i,17) = num2cell(smallwtfrac);
    noseps_props(i,34)= num2cell(largewtfrac);
end

toc
end