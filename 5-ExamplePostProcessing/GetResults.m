clear; clc; close all;

FT = 'F:\Stephen\2019-Torsion\2019-2-1 Torsion Full';
FB = 'F:\Stephen\2019-Bending\2019-2-12 Bending Full';
FU = 'F:\Stephen\2019-Bending\2019-2-12 Bending Full\Uniaxial Tension Tests';

txt_list_T = dir([FT,'\*.txt']);
txt_list_B = dir([FB,'\*.txt']);
txt_list_U = dir([FU,'\*.txt']);
txt_list=txt_list_B;

%Preallocate Results;
[seed,L,R,Eb,Et,Eu,Ot,Ob]=deal(NaN(length(txt_list),1));
for i = 1:length(txt_list)
   % Setup File Paths
   fileT = [FT '\' txt_list(i).name];
   fileB = [FB '\' txt_list(i).name];
   fileU = [FU '\' txt_list(i).name];
   fileM = [fileB(1:end-4) '.mat'];
   % Load In Mat File
   load(fileM,'final_nodes','final_fibers','l_fiber','index_inp');
   seed(i) = index_inp;
   L(i) = round(range(final_nodes(:,2)));
   R(i) = round(range(final_nodes(:,3)))/2;
   I = (R(i))^4*pi/4;
   J = (R(i))^4*pi/2;
   % === === Torsion === ===
       RawData = importdata(fileT);
       FullData=array2table(RawData.data,'VariableNames',RawData.colheaders);
       POS_i = [FullData.X FullData.Y FullData.Z];
       POS_f = POS_i + [FullData.U1 FullData.U2 FullData.U3];
       u=POS_i(1,:); %Vectors must be 3D
       v=POS_f(1,:);
       theta = atan2(norm(cross(u,v)),dot(u,v));
       % Calculate Applied Moment
       M_arms = POS_f; %Find Moment Arms
       Moment=0;
       for j=1:height(FullData)
          Moment=Moment+M_arms(j,3)*FullData.RF2(j);
          Moment=Moment-M_arms(j,2)*FullData.RF3(j);
       end
       % Calculate Modulus
       G = Moment*L(i)/(2*theta)./J; %Note this is Shear Modulus
       Et(i) = G * 2 * (1+0.33); % assumed poisson
   % === === Bending === ===
       RawData = importdata(fileB);
       FullData=array2table(RawData.data,'VariableNames',RawData.colheaders);
       POS_i = [FullData.X FullData.Y];
       POS_f = POS_i + [FullData.U1 FullData.U2];
       % Find Center before and after load
       tform = fitgeotrans(POS_i,POS_f,'nonreflectivesimilarity');
       Center_i = [L(i)/2,0];
       Center_f = transformPointsForward(tform,Center_i);
       % Find Angle Change
       Tip_i = [L(i)/2,R(i)];
       Tip_f = transformPointsForward(tform,Tip_i);
       u=[Tip_i-Center_i,0]; %Vectors must be 3D
       v=[Tip_f-Center_f,0];
       theta = atan2(norm(cross(u,v)),dot(u,v));
       % Calculate Applied Moment
       M_arms = (POS_f-Center_f); %Find Moment Arms
       Moment=0;
       for j=1:height(FullData)
          Moment=Moment+M_arms(j,1)*FullData.RF2(j);
          Moment=Moment-M_arms(j,2)*FullData.RF1(j);
       end       
       % Calculate Modulus
       rho = L(i)./2./theta;
       %https://www.engineeringtoolbox.com/cantilever-beams-d_1848.html
       Eb(i)  = - Moment.*rho./I; %Bending Modulus
   % === === Tension === ===
       RawData = importdata(fileU);
       FullData=array2table(RawData.data,'VariableNames',RawData.colheaders);
       POS_i = [FullData.X FullData.Y];
       POS_f = POS_i + [FullData.U1 FullData.U2];
       RM_vec = [FullData.RF1,FullData.RF2];
       A = pi*R(i)^2; %Assumed Cylinder
       % Calculate Values
       F = sum(FullData.RF1);
       ex = mean(FullData.U1)/(L(i)/2);
       Eu(i) = F/A/ex;  %Modulus from uniaxial tension
   % === === Save Omega Vals === ===
   Ot(i) = Et(i)/Eu(i);
   Ob(i) = Eb(i)/Eu(i);
end

data=array2table([seed,L,R,Eb,Et,Eu,Ot,Ob],'VariableNames',{'seed','L','R','Eb','Et','Eu','Ot','Ob'});
save('FullScaleResults','data');