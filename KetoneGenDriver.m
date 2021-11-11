% Edited by Jake Miller, 4/20/2021. Complete species list (in order):
% 1) Acetic acid
% 2) Propanoic acid
% 3) 2-methyl propanoic acid
% 4) Butanoic acid
% 5) (2?)-methyl butanoic acid
% 6) Pentanoic acid
% 7) Hexanoic acid
% 8) Heptanoic acid
% 9) Octanoic acid
% 10) isopropanol
% 11)2-butanol
% 12)3-methyl-2-butanol
% 13)2-pentanol
% 14)4-methyl-2-pentanol
% 15)2-hexanol
% 16)2-heptanol
% 17)2-octanol
% 18)2-nonanol
% 19)3-pentanol
% 20)2-methyl-3-pentanol
% 21)3-hexanol
% 22)5-methyl-3-hexanol
% 23)3-heptanol
% 24)3-octanol
% 25)3-nonanol
% 26)3-decanol
% 27)2,4-dimethyl 3-pentanol
% 28)2-methyl 3-hexanol
% 29)2,5-dimethyl 3-hexanol
% 30)2-methyl 3-heptanol
% 31)2-methyl 3-octanol
% 32)2-methyl 3-nonanol
% 33)2-methyl 3-decanol
% 34)4-heptanol
% 35)2-methyl 4-heptanol
% 36)4-octanol
% 37)4-nonanol
% 38)4-decanol
% 39)4-undecanol
% 40)2,6-dimethyl 4-heptanol
% 41)2-methyl 4-octanol
% 42)2-methyl 4-nonanol
% 43)2-methyl 4-decanol
% 44)2-methyl 4-undecanol
% 45)5-nonanol
% 46)5-decanol
% 47)5-undecanol
% 48)5-dodecanol
% 49)6-undecanol
% 50)6-dodecanol
% 51)6-tridecanol
% 52)7-tridecanol
% 53)7-tetradecanol
% 54)8-pentadecanol

function [Cketone, ketoneend,light, heavy] = KetoneGenDriver(VFA_molfracs)
global k n total Nacid Nketone
% Acid reactant carbon number
Cacid = [2 3 4 4 5 5 6 7 8]';
n=length(Cacid);
% Ketone product carbon numbers
for i =1:n
        for j=i:n
        l=Cacid(i)+Cacid(j)-1;
        Cketone(l)=l;
        end         
end
ketone0= 2*Cacid(1)-1;
Cketone=Cketone(:,ketone0:end)';
Sketone=length(Cketone);
%Position of iso acid:
m=1;
for i =1:n-1
    if Cacid(i)== Cacid(i+1)
        iso(m)=i;
        m=m+1;
    end        
end


Nacid = n; % number of acid reactants
Nketone = n*(n+1)/2; %number of ketone products
total = Nacid +Nketone; % total number of acids and ketones
C0 = zeros(total,3);
%Comp_list = {'Comp1','Comp2', 'Comp3', 'Comp4', 'Comp5', 'Comp6', 'Comp7','Comp8','Comp9','Comp10','Comp11'};
%run = Comp_list{3};
%data = [1 0 0 0 0 0 0 0 0 0 0];
%[ratio,data] = KetoneGenData(run);
% H = ratio(1);
C0 = VFA_molfracs; %initial acid and ketone concentration

for i = 1:n
    C0(i)=VFA_molfracs(i);
end
for i = n+1:total
    C0(i)= 0;
end
%rate constants
k = zeros (n,n);
  for i = 1:n
      for j = i: n
    k (j,i) = 0.005; 
    if i==j
        k(j,i)= k(j,i)/2;
    end
      end
  end
%Solve ODES  
tspan = [0 1000000];
% options= odeset('RelTol',1E-10,'AbsTol',1E-10);
[t,C]=ode45(@ODEs_KetoneGen,tspan,C0); % 77 row of t
% acid profile
acid = C(:,1:n);
% ketone profile
ketone = C(:,n+1:total);
tmax = length(t);


% Sketone = 2*n-1;
% Cketone = zeros (Sketone,1);
% Cketone(1) = 2*Cacid (1)-1;
% for i = 2:Sketone
%     Cketone(i) = Cketone(i-1) +1;
% end
% final ketone concentration vector
ketonefinal = ketone(tmax,:)';
ketoneend=ketone(tmax,:)';
% Create 2D array out of the ketonefinal vector
for i = 1:n-1 
z= zeros (i,1);
ketonefinal = [ketonefinal(1:n*i);z; ketonefinal(n*i+1:end)]; % Add more rows to ketonefinal to prepare to reshape 1D to 2D array
end
ketone2D = reshape(ketonefinal',n,n); % Create 2D array out of 1D array
% Group ketones by carbon numbers
Ketone_branch=zeros(n,15);
Ketone_normal=zeros(n,15);
for i=1:n
    if ismember(i,iso)==1
        for j=i:n
           Cketone_branch = Cacid(i)+Cacid(j)-1;
           Ketone_branch(i,Cketone_branch)=ketone2D(j,i); 
        end
    else
    for j=i:n
        if ismember(j,iso)==1
            Cketone_branch = Cacid(i)+Cacid(j)-1;
            Ketone_branch(i,Cketone_branch)=ketone2D(j,i);
        else                       
        Cketone_linear = Cacid(i)+Cacid(j)-1;
        Ketone_linear(i,Cketone_linear)=ketone2D(j,i);
        end
        
    end      
        
    end
end
Branch_Ketone=sum(Ketone_branch(:,Cketone))';
Linear_Ketone =sum(Ketone_linear(:,Cketone))';
Ketone = Branch_Ketone+Linear_Ketone;
light = [];
heavy = [];
%Create a vector of <C8 ketone or alcohols
for i =1:n
    for j =i:n
        Calcohol = Cacid(i)+Cacid(j)-1;
        if Calcohol < 8
            light =[light; ketone2D(j,i)];
        else
            heavy =[heavy;ketone2D(j,i)];
        end
    end
end


% for i = 1:l
%     if light(i)==0
%         light(i)=[];
%     end
% end
% for i = 1:h
%     if heavy(i)==0
%         heavy(i)=[];
%     end
% end
%         
% light
% heavy
% ketone2D;
        
% ketone2D_new = rot90(ketone2D);
% % Group ketones by carbon numbers
%  for i = 1:Sketone
%     Ketone(Sketone+1-i,:) = sum (diag(ketone2D_new,n-i));
%  end
% Molecular weight calculation
MW = Cketone.*12.0107+ (Cketone.*2+2)*1.00784;
MW_acid = Cacid.*12.0107 + Cacid.*2*1.00784 + 2*16;
MW_ketone = Cketone.*12.0107 + Cketone.*2*1.00784 +16;
%From Glenn experiment
% MassKetoneExp = [0
% 0.91
% 2.78
% 3.70
% 11.34
% 7.11
% 26.05
% 5.57
% 28.03
% 0.00
% 10.29
% ];
% MoleKetoneExp = MassKetoneExp./MW_ketone;
% Ketone =MoleKetoneExp
% Mole and mass fraction of Ketone
% MassHC= Ketone.*MW;
% 
% for i = 1:Sketone
%     MoleFracKetone(i) = Ketone(i)/sum(Ketone)*100;
%     MassFracKetone (i) = Ketone (i)*MW(i)/sum(MassHC)*100;
% end
for i = 1:Sketone
    MoleFracKetone_Linear(i) = Linear_Ketone(i)/sum(Ketone)*100;
    MoleFracKetone_Branch(i) = Branch_Ketone(i)/sum(Ketone)*100;
    MoleFracKetone(i) = Ketone(i)/sum(Ketone)*100;
   end
MoleFracKetone_final= [MoleFracKetone_Linear; MoleFracKetone_Branch]';
% %%% Boiling Point Profile
% % Boiling point critera
% CummJet =[5 10 30 50 70 90 100];
% JetMin = [100 130 160 180 200 225 250];
% JetMax = [140 160 185 210 230 260 300];
% % Boiling point of hydrocarbon products
% 
% % MW = [114.23 128.2 142.29 156.31 170.33 184.37 198.39 212.42];
% NormBP = [126 151 174 196 215 234 254 270];
% NormIsoBP = NormBP-7;
% 
% for i = 1:Sketone
% if Cketone(i)==8
% SelectCHC = Cketone (i:Sketone);
% SelectHC = Ketone(i:Sketone);
% SelectMW = MW(i:Sketone);
% break
% end
% end
% SelectMassHC =SelectHC.*SelectMW;
% TotalMass = sum(SelectMassHC);
% MassFracHC = SelectMassHC/TotalMass*100;
% MassFracHC_n = MassFracHC*0.85;
% MassFracHC_iso = MassFracHC*0.15;
% MassFracHC_final = [MassFracHC_n, MassFracHC_iso];
% 
% density = 730; % unit kg/m^3
% % Cummulative Mass of HC by carbon number
% for i = 1:size(MassFracHC)
%     if i ==1
%     CummMass(i) = MassFracHC(i);
%     else
%         CummMass (i) = MassFracHC(i)+ CummMass (i-1);
%     end
% end
% 
% % Cummulative Mass of HC by carbon number and iso-normal paraffins
% MassFracHC_new = [MassFracHC_iso.';MassFracHC_n.'];
% MassFracHC_new = MassFracHC_new(:);
% SelectCHC_new = [SelectCHC.';SelectCHC.'];
% SelectCHC_new = SelectCHC_new(:);
% for i = 1:size(MassFracHC_new)
%     if i ==1
%     CummMass_new(i) = MassFracHC_new(i);
%     else
%         CummMass_new (i) = MassFracHC_new(i)+ CummMass_new (i-1);
%     end
% end
% 
% % Carbon distribution to jet, naphtha
% Carbon_HC = Cketone.*Ketone
% for i = 1:Sketone
% if Cketone(i) <8
%     NapCarbon = sum(Carbon_HC(1:i,:))
%     nap_mass = Ketone(1:i,:).*MW(1:i,:);
%     NapMass = sum(nap_mass);
%     else
%     JetCarbon = sum(Carbon_HC(i:end,:));
%     jet_mass = Ketone(i:end,:).*MW(i:end,:);
%     JetMass = sum (jet_mass);
%     break
% end
% end
% % For Glenn experimnet
% Co2Carbon = sum(Ketone);
% TotalCarbonInitial = NapCarbon + JetCarbon + Co2Carbon;
% TotalMassInitial = NapMass + JetMass + Co2Carbon *44.01 + Co2Carbon*18;
% % Co2Carbon = sum(ketonefinal);
% NapCarbon
% JetCarbon
% % AcidCarbonFinal=sum(Cacid.*acid(tmax,:)')
% % TotalCarbonInitial = sum(Cacid.*data.')
% % TotalCarbonFinal = NapCarbon +JetCarbon +Co2Carbon + AcidCarbonFinal
% NapCarbonMolFrac = NapCarbon/TotalCarbonInitial*100;
% JetCarbonMolFrac = JetCarbon/TotalCarbonInitial*100;
% Co2CarbonMolFrac = Co2Carbon/TotalCarbonInitial*100;
% %Mass distribution to jet, naptha, CO2
% % TotalMassInitial = sum(MW_acid.*data.');
% NapMassFrac = NapMass/TotalMassInitial*100;
% JetMassFrac = JetMass/TotalMassInitial*100;
% Co2MassFrac = Co2Carbon*44.01/TotalMassInitial*100;

% Visualization

% % Boiling point plots
% plot (CummMass,NormIsoBP,'-o')
% plot(CummJet,JetMin,'-o')
% plot(CummJet,JetMax,'-o')
% Distribution of ketone products
%%%Bar graph to show mole distribution by C number 
% figure
% hold on
% bar(Cketone,MoleFracKetone)
% labels = arrayfun(@(value) num2str(value,2),MoleFracKetone,'UniformOutput',false);
% text(Cketone,MoleFracKetone,labels,...
%   'HorizontalAlignment','center',...
%   'VerticalAlignment','bottom')
% ylabel('Mole Percentage')
% xlabel('Product Carbon Number')

%%%Bar graph to show mass distribution by C number and linear, branch ketones
% figure
% hold on
% bar(Cketone,MoleFracKetone_final,'stacked')
% labels = arrayfun(@(value) num2str(value,2),MoleFracKetone,'UniformOutput',false);
% legend ('Linear','Branch');
% text(Cketone,MoleFracKetone,labels,...
%   'HorizontalAlignment','center',...
%   'VerticalAlignment','bottom')
% ylabel('Mole Percentage')
% xlabel('Ketone Carbon Number')

% %%%Bar graph to show mass distribution by C number and normal, iso paraffins
% figure
% hold on
% bar(SelectCHC,MassFracHC_final,'stacked')
% labels = arrayfun(@(value) num2str(value,2),MassFracHC,'UniformOutput',false);
% legend ('n-paraffins','isoparaffins');
% text(SelectCHC,MassFracHC,labels,...
%   'HorizontalAlignment','center',...
%   'VerticalAlignment','bottom')
% ylabel('Mass Percentage')
% xlabel('Paraffin Carbon Number')
% 
% %%% Plot to show cummulative mass by carbon number and iso/normal
% figure
% hold on
% x= linspace(1,16,16);
% plot(x,CummMass_new,'-*')
% 
% grid on
% ylabel('Cummulative Mass (%)')
% xlabel('Paraffin Product')
% xticks(1:1:24)
% xticklabels({'iso-C8','n-C8','iso-C9','n-C9','iso-C10','n-C10','iso-C11','n-C11','iso-C12','n-C12','iso-C13','n-C13','iso-C14','n-C14','iso-C15','n-C15',})
% xtickangle(45)

% cl='rgbmc'
% for j=1:numel(MoleFracKetone)
%   bar(j,MoleFracKetone(j),cl(MoleFracKetone(j)))
%   hold on
% end
% Time profile of ketone products
% hold on
% subplot(2,1,1)
% hold on
% plot(t,ketone(:,1),'o')
% plot(t,ketone(:,6),'*')
% plot(t,ketone(:,10),'^')
% plot(t,ketone(:,13),'+')
% plot(t,ketone(:,15),'-')
% subplot (2,1,2)
% plot(t,C(:,20),'+')
% plot(t,C(:,9),'^')
% plot(data(1,:),data(2,:),'*m')
% plot(data(1,:),data(3,:),'*b')
% plot(data(1,:),data(4,:),'*r')
% plot(data(1,:),data(5,:),'*k')
% legend('c_F_A','c_E_t_h_e_r','c_L_A','c_B_P')
% xlabel('Reaction Time [hr]')
% ylabel('Concentration [mol/L]')



end
