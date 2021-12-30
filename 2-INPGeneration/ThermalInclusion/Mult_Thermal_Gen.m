clear; clc; close all;

% Generate Folder
Date = char(datetime('now','Format','yyyy-MM-dd'));
SimTitle = 'Rigid Inclusion Displacement Tests';
Dir = [Date ' ' SimTitle];
clear Date SimTitle
mkdir(Dir)

% User Inputs
Net_Fold = 'C:\Users\Stephen\Desktop\NetGenOpt\V2';
Net_Base = 'ContractionDomain_500x500x2_Seed_';
Net_Seed = 1001;

Radius_domain = 250;
Radius_particle = 15;
Radius_inclusion = Radius_particle;
Percent_Contract = 20; 

Inc_Distance = 50:50:200;
Inc_Angle = 0:pi/2:3/2*pi;

% INP Generation
for iSeed=Net_Seed
   Net_Path = [fullfile(Net_Fold,[Net_Base,num2str(iSeed)],[Net_Base,num2str(iSeed)]),'.mat'];
   File_Path = Cut_2_Annulus(Dir,Net_Path,Radius_particle,Radius_domain);
   for iAng = Inc_Angle
   for iGap = Inc_Distance
   for iStrain = Percent_Contract
       ThermalShrink_INP(File_Path,Radius_particle,Radius_inclusion,iSeed,iAng,iGap,iStrain)
   end
   end
   end   
end
