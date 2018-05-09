%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is part of StabFem Project.
%   (set of programs in Freefem to perform global stability calculations  
%   in hydrodynamic stability)
%
%	File: Main_VIV_Free_Movement.m
%   Contributours: David Fabre, Diogo Sabino ...
%   Last Modification: Diogo Sabino, 27 April 2018
%   
%	This script models a flow aroung a rigid circular cilinder mounted
%	in a spring-damped system to simulate VIV. Case
%
%   To give a good use of this script, run each section at a time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh creation and convergence using adapt mesh
%clear all
global ffdataharmonicdir verbosity

close all
%baseflow.mesh.problemtype='2D'
run('../SOURCES_MATLAB/SF_Start.m');

disp(' ');  disp(' GENERATING  MESH : [-40:120]x[0:40] '); disp(' ');


baseflow=SF_Init('Mesh_Cylinder.edp', [-40 80 40]);
baseflow=SF_BaseFlow(baseflow,'Re',1);
baseflow=SF_BaseFlow(baseflow,'Re',10);
baseflow=SF_BaseFlow(baseflow,'Re',60);

fig_mesh=plotFF(baseflow,'mesh'); %pause;
set(fig_mesh, 'HandleVisibility', 'off'); %IMFT just have 3 licences, so...

disp(' ');disp('ADAPTING MESH FOR RE=60 ');disp(' ');
baseflow=SF_Adapt(baseflow,'Hmax',10,'InterpError',0.005);

plotFF(baseflow,'mesh');%pause(0.1);
baseflow=SF_BaseFlow(baseflow,'Re',60);

% 
disp(' ');disp('ADAPTING MESH FOR RE=60 ACORDING TO EIGENVALUE ');disp(' ');
% adaptation du maillage sur un mode propre
[ev,em] = SF_Stability(baseflow,'shift',0.04+0.74i,'nev',1,'type','D');
[baseflow,em]=SF_Adapt(baseflow,em,'Hmax',10,'InterpError',0.02); 

plotFF(baseflow,'mesh');%pause(0.1);

%% PARTIE VIV Search in the spectral space

%disp('Press any key to continue to VIV Case');
%pause;
%To pass to the eigenvalue problem taking into acount the VIV
%baseflow.mesh.problemtype='2D_VIV';

%% Validation Phase: Parameter Definition
Re=60; verbosity=1;
baseflow = SF_BaseFlow(baseflow,'Re',Re); %Do it with problemtype 2D
baseflow.mesh.problemtype='2D_VIV'; verbosity=10;

m_star=5;
mass=m_star*pi/4;
%mass=m_star;
U_star=[3.1:0.1:8.5 8.5:0.05:11 ];%Pour le mode CYL%
%Pour le mode FLUID%U_star=[4:0.3:7 7:0.05:11 ];
%U_star=[4:0.3:20];%Pour le mode FLUID 2nd essay%


Stiffness_table=pi^3*m_star./((U_star).^2) ;
%Stiffness_table = [1.5:0.5:10 10: 1 :20];
%Stiffness_table = [10: 10 :100];
sigma_tab = []; verbosity=10;

%% Validation Phase: See spectrum
for STIFFNESS=Stiffness_table
    %STIFFNESS=Stiffness_table(8) %for searching just for one k
    realshift=-0.1;
    imshift=1.8:-0.2:1.8;
    for ii=realshift
        for jj=imshift
            [ev,em] = SF_Stability(baseflow,'shift',ii+jj*1i,'nev',5,'type','D','STIFFNESS',STIFFNESS,'MASS',mass,'DAMPING',0,'Frame','A','PlotSpectrum','yes');
        end
    end
end

%% Validation Phase: Follow a mode along the Stiffness_table to reproduze Navrose's images

%Manually put the right shift
%Pour le MODE CYL mstar20
%[ev,em] = SF_Stability(baseflow,'shift',-0.03+1.5i,'nev',1,'type','D','STIFFNESS',Stiffness_table(1),'MASS',mass,'DAMPING',0,'Frame','A','PlotSpectrum','yes');
%Pour mode fluid mstar20
sigma_tab = [];
[ev,em] = SF_Stability(baseflow,'shift',-0.11+1.81i,'nev',1,'type','D','STIFFNESS',Stiffness_table(1),'MASS',mass,'DAMPING',0,'Frame','A','PlotSpectrum','yes');


for STIFFNESS=Stiffness_table
	disp('Next sigma calculated'); disp(STIFFNESS);
	[ev,em] = SF_Stability(baseflow,'shift','cont','nev',1,'type','D','STIFFNESS',STIFFNESS,'MASS',mass,'DAMPING',0,'Frame','A','PlotSpectrum','yes'); 
    sigma_tab = [sigma_tab ev];
end


 figure;hold on;
 plot(U_star,real(sigma_tab),'r*');
 title('amplification rate');
% %
 figure;hold on;
 plot(U_star,imag(sigma_tab),'b*');
 title('oscillation rate');
 
%save('modeCYL_zoom.mat','sigma_tab','Stiffness_table','U_star');




%% Reynolds Analysis
%Spectrum exploration
Re_tab = 30: 10 : 60;
for Re=Re_tab
    verbosity=1;
    baseflow = SF_BaseFlow(baseflow,'Re',Re);
    verbosity=10;
    realshift=0;
    imshift=0.1:0.1:1;
    
    for ii=realshift
        for jj=imshift
            [ev,em] = SF_Stability(baseflow,'shift',ii+jj*1i,'nev',10,'type','D','STIFFNESS',0.1,'MASS',1,'DAMPING',0.01,'Frame','A','PlotSpectrum','yes');
        end
    end
end

%Stability Branch
baseflow=SF_BaseFlow(baseflow,'Re',60);
[ev,em] = SF_Stability(baseflow,'shift',0.04678+0.7451i,'nev',1,'type','D','STIFFNESS',0.1,'MASS',1,'DAMPING',0.01,'Frame','A','PlotSpectrum','yes');
Re_tab = 60: -2.5 : 30;
sigma_tab = [];

for Re = Re_tab
	verbosity=1;
    baseflow = SF_BaseFlow(baseflow,'Re',Re);
    verbosity=10;
	[ev,em] = SF_Stability(baseflow,'shift','cont','nev',1,'type','D','STIFFNESS',0.1,'MASS',1,'DAMPING',0.01,'Frame','A','PlotSpectrum','yes'); 
    sigma_tab = [sigma_tab ev];
end

 figure;hold on;
 plot(Re_tab,real(sigma_tab),'r*');
 title('amplification rate');
% %
 figure;hold on;
 plot(Re_tab,imag(sigma_tab),'b*');
 title('oscillation rate');


%boucle des tests de l'apres-midi du jour D21: 3/05/2018
if(1==0)
    baseflow.mesh.problemtype='2D_VIV';
    Re=60;
    baseflow = SF_BaseFlow(baseflow,'Re',Re);
    verbosity=10;
    
    realshift=0.07:-0.01:-0.03;
    %realshift=0.02;
    %imshift=0.76:-0.01:0.72;
    imshift=0.76;
    
    for ii=realshift
        for jj=imshift
            [ev,em] = SF_Stability(baseflow,'shift',ii+jj*1i,'nev',5,'type','D','STIFFNESS',77.5157,'MASS',15.7080,'DAMPING',0,'Frame','R','PlotSpectrum','yes');
        end   
    end   
end







    
%% starting point ; autres choses 
if(1==0)

[ev,em] = SF_Stability(baseflow,'shift',-.03+.72i,'nev',1,'type','D','STIFFNESS',1,'MASS',30,'DAMPING',0,'Frame','A');
 



for Re = Re_tab
    verbosity=1;
    baseflow = SF_BaseFlow(baseflow,'Re',Re);
    verbosity=10;
    [ev,em] = SF_Stability(baseflow,'shift',-.03+.72i,'nev',1,'type','D','STIFFNESS',0.1,'MASS',1,'DAMPING',0.01,'Frame','A');
    
    disp('___________________________________________________')%To separete iters
end

figure;hold on;
plot(Re_tab,real(sigma_tab),'r*');
title('amplification rate');
%
figure;hold on;
plot(Re_tab,imag(sigma_tab),'b*');
title('oscillation rate');


%% autres choses qui avait dans le document cas forcé 
% 
% 
% 
% end
% if(1==0)
% 
% disp(' ');
% disp('COMPUTING STABILITY BRANCH FOR VIV CASE (M=10,K=1) ')
% disp(' ');
% % starting point
% baseflow=SF_BaseFlow(baseflow,'Re',40);
% [ev,em] = SF_Stability(baseflow,'shift',-0.015851+0.75776i,'nev',1,'type','D','STIFFNESS',1,'MASS',30,'DAMPING',0);
% 
% Re_tab = [40 : 2.5 : 80];
% sigma_tab = [];
% for Re = Re_tab
%     baseflow=SF_BaseFlow(baseflow,'Re',Re);
%     [ev,em] = SF_Stability(baseflow,'shift','cont','nev',1,'type','D','STIFFNESS',1,'MASS',10,'DAMPING',0);
%     sigma_tab = [sigma_tab ev];
% end
% 
% figure(1); hold on;
% plot(Re_tab,real(sigma_tab),'r--');
% title('amplification rate');
% 
% figure(2);hold on;
% plot(Re_tab,imag(sigma_tab),'b--');


%% Choses que David avait fait dans le script VIV test
% Matlab interface for GlobalFem
%
%   (set of programs in Freefem to perform global stability calculations in hydrodynamic stability)
%  
%  THIS SCRIPT DEMONSTRATES THE MESH CONVERGENCE USING ADAPTMESH 
%   TO BASEFLOW AND EIGENMODE.

global ffdataharmonicdir verbosity

close all

run('../SOURCES_MATLAB/SF_Start.m');

disp(' ');  disp(' GENERATING  MESH : [-40:120]x[0:40] '); disp(' ');

baseflow=SF_Init('Mesh_Cylinder.edp', [-40 80 40]);
baseflow=SF_BaseFlow(baseflow,'Re',1);
baseflow=SF_BaseFlow(baseflow,'Re',10);
baseflow=SF_BaseFlow(baseflow,'Re',60);
fig_mesh=plotFF(baseflow,'mesh'); %pause;
set(fig_mesh, 'HandleVisibility', 'off'); %IMFT just have 3 licences, so...

disp(' ');disp('ADAPTING MESH FOR RE=60 ');disp(' ');
baseflow=SF_Adapt(baseflow,'Hmax',10,'InterpError',0.005);
plotFF(baseflow,'mesh');%pause(0.1);
baseflow=SF_BaseFlow(baseflow,'Re',60);
% 
disp(' ');disp('ADAPTING MESH FOR RE=60 ACORDING TO EIGENVALUE ');disp(' ');
% adaptation du maillage sur un mode propre
[ev,em] = SF_Stability(baseflow,'shift',0.04+0.74i,'nev',1,'type','D');
[baseflow,em]=SF_Adapt(baseflow,em,'Hmax',10,'InterpError',0.005); %j'ai
%changé ca!!!!
plotFF(baseflow,'mesh');%pause(0.1);
% 
% what is this ?
% %disp(' ');
% %disp(' THRESHOLD : ')
% %[baseflowC,emC]=SF_FindThreshold(baseflow,em);



% VIV TEST
close all;

disp(' ');
disp('COMPUTING STABILITY BRANCH (fixed cylinder) ')
disp(' ');
% starting point
baseflow=SF_BaseFlow(baseflow,'Re',40);
[ev,em] = SF_Stability(baseflow,'shift',-0.03+.72i,'nev',1,'type','D');



Re_tab = [40 : 5 : 80];
sigma_tab = [];
for Re = Re_tab
    baseflow = SF_BaseFlow(baseflow,'Re',Re);
    [ev,em] = SF_Stability(baseflow,'shift','cont','nev',1,'type','D');
    sigma_tab = [sigma_tab ev];
end
figure(1);hold on;
plot(Re_tab,real(sigma_tab),'r*');
title('amplification rate');

figure(2);hold on;
plot(Re_tab,imag(sigma_tab),'b*');
title('oscillation rate');
pause(0.1);




%END. old scripts new



disp(' ');
disp('COMPUTING STABILITY BRANCH FOR VIV CASE (M=10,K=1) ')
disp(' ');
% starting point
baseflow=SF_BaseFlow(baseflow,'Re',40);
[ev,em] = SF_Stability(baseflow,'shift',-0.015851+0.75776i,'nev',1,'type','D','STIFFNESS',1,'MASS',30,'DAMPING',0);

Re_tab = [40 : 2.5 : 80];
sigma_tab = [];
for Re = Re_tab
    baseflow=SF_BaseFlow(baseflow,'Re',Re);
    [ev,em] = SF_Stability(baseflow,'shift','cont','nev',1,'type','D','STIFFNESS',1,'MASS',10,'DAMPING',0);
    sigma_tab = [sigma_tab ev];
end

figure(1); hold on;
plot(Re_tab,real(sigma_tab),'r--');
title('amplification rate');

figure(2);hold on;
plot(Re_tab,imag(sigma_tab),'b--');

end