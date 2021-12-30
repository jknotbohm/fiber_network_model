clear; clc; close all;

%Provide path to the tension results .txt file
%this should come from the matlab ODB_NODE file
fileU= uigetfile;

% === === Tension === ===
   RawData = importdata(fileU);
   FullData=array2table(RawData.data,'VariableNames',RawData.colheaders);
   POS_i = [FullData.X FullData.Y];
   POS_f = POS_i + [FullData.U1 FullData.U2];
   RM_vec = [FullData.RF1,FullData.RF2];
   R=max(FullData.Y);
   L=2*max(FullData.x);
   A = pi*R^2; %Assumed Cylinder
   % Calculate Values
   F = sum(FullData.RF1);
   ex = mean(FullData.U1)/(L/2);
   E = F/A/ex;  %Modulus from uniaxial tension