%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is part of StabFem Project.
%   (set of programs in Freefem to perform global stability calculations  
%   in hydrodynamic stability)
%
%	File: Main_VIV_Free_Movement.m
%   Contributours: David Fabre, Diogo Sabino ...
%   Last Modification: Diogo Sabino, 16 May 2018
%   
%	This script models a flow aroung a rigid circular cilinder mounted
%	in a spring-damped system to simulate VIV.
%
%   To give a good use of this script, run each section at a time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh: creation and convergence using adapt mesh
%clear all
close all
global ffdataharmonicdir verbosity

run('../SOURCES_MATLAB/SF_Start.m');
%CHOOSE the domain parameters:
domain_parameters=[-50 50 50];
disp_info=[' GENERATING  MESH :[' num2str(domain_parameters(1)) ':' num2str(domain_parameters(2)) ']x[0:' num2str(domain_parameters(3)) ']'];
disp(' ');  disp(disp_info); disp(' ');
verbosity=10;
baseflow=SF_Init('Mesh_Cylinder.edp', domain_parameters);
baseflow=SF_BaseFlow(baseflow,'Re',1);
baseflow=SF_BaseFlow(baseflow,'Re',10);
baseflow=SF_BaseFlow(baseflow,'Re',60);

%fig_mesh=plotFF(baseflow,'mesh'); %pause;
%set(fig_mesh, 'HandleVisibility', 'off'); %IMFT just have 3 licences, so...

disp(' ');disp('ADAPTING MESH FOR RE=60 ');disp(' ');
baseflow=SF_Adapt(baseflow,'Hmax',10,'InterpError',0.005);

%plotFF(baseflow,'mesh');%pause(0.1);
baseflow=SF_BaseFlow(baseflow,'Re',60); % not needed, I think

disp(' ');disp('ADAPTING MESH FOR RE=60 ACORDING TO EIGENVALUE ');disp(' ');
% adaptation du maillage sur un mode propre
baseflow.problemtype='2D'; %Used just to save time, using Stab_2D.edp
[ev,em] = SF_Stability(baseflow,'shift',0.04+0.74i,'nev',1,'type','D');
baseflow.problemtype='2D_VIV';
[baseflow,em]=SF_Adapt(baseflow,em,'Hmax',10,'InterpError',0.02);
%plotFF(baseflow,'mesh');%pause(0.1);
%[baseflow,em]=SF_Adapt(baseflow,em,'Hmax',5,'InterpError',0.01);
%plotFF(baseflow,'mesh');%pause(0.1);
%[baseflow,em]=SF_Adapt(baseflow,em,'Hmax',5,'InterpError',0.01); 
%plotFF(baseflow,'mesh');%pause(0.1);
%baseflow=SF_Split(baseflow);
plotFF(baseflow,'mesh');%pause(0.1);

%CHOOSE folder for saving data:
General_data_dir='./Final_Results_v20/'; %a array of char

domain_identity={[ num2str(domain_parameters(1)) '_' num2str(domain_parameters(2)) '_' num2str(domain_parameters(3)) '/']};
%CHOOSE the name of the folder mesh: (It's good to choose the name because we can adapt mesh during)
mesh_identity={'Adapt_mode_Hmax10_InterError_0.02/'};
savedata_dir={[ General_data_dir domain_identity{1} mesh_identity{1}]};
%savedata_dir=[General_data_folder domain_identity mesh_identity];

%% Parameters' Definition
%CHOOSE the Re to test:
Re=23; verbosity=10;
baseflow = SF_BaseFlow(baseflow,'Re',Re);

%CHOOSE the m_star to test:
m_star_tab=[10 20 40 100 1000]
%m_star=10;
for m_star=m_star_tab
    mass=m_star*pi/4;
    %U_star=[11:-0.2:9.5 9.5:-0.05:6.5 6.5:-0.2:3];
    U_star=[3:0.2:6.5 6.5:0.05:11 ];%11:0.1:15
    %U_star=[3:0.2:4.5];% TOTO efacer
    
    Stiffness_table=pi^3*m_star./((U_star).^2);
    
    % Spectrum search: See spectrum if wanted, for discover the shift
    if(1==0)
        sigma_tab=[];
        %Re and mass are already defined above. Define numbers or tables for:
        STIFFNESS_to_search=Stiffness_table(1); %Use whatever we want
        RealShift=[0.05 0 -0.03];
        ImagShift=[0.75];
        
        nev=5; %Normally here we don't use just one
        [baseflow,sigma_tab] = SF_FreeMovement_Spectrum('search',baseflow,sigma_tab,RealShift,ImagShift,STIFFNESS_to_search,mass,nev);
        filename={'01spectrum_search'};
        SF_Save_Data('spectrum',General_data_dir,savedata_dir,Re,m_star,filename,Stiffness_table,U_star,sigma_tab);
    end
    % Mode Follow: Follow a mode along the Stiffness_table/U_star
    
    %Manually put the right shift, acordding to the above spectrum(
    %For the pair (Re;[m_star])(starting at U_star(1)=3 )
    %For the STRUCTURE Mode:
    %shift=0+2.1i:                                                                      %%50;[100])
    %shift=0+2i:	([60];[20,10])   ([40-19];[1000-10])                       %%(15;[100])
    %shift=0+1.8i:  (60;[5])              (40;[4.73])
    %shift=0+1i:                                                                         %%(40;[0.7])
    
    %For the FLUID Mode :
    %shift=0.05+0.75i:	(60;[20,10,5])
    %shift=-0.03+0.75i:               (40,[300-5])
    %shift=-0.07+0.62i:                                                                 %%%%(21;[10])
    %close all
    
    %CHOOSE shift:
    RealShift=0; ImagShift=2; %Normally, for Structure
    %RealShift=-0.03; ImagShift=0.75; %Normally, for Fluid
    
    %CHOOSE the one for save data w/ a good name:
    modename={'02modeSTRUCTURE'};
    %modename={'03modeFLUID'};
    
    
    sigma_tab=[];
    nev=1; %Normally if's just one, but if shift is wrong, it helps put more
    
    [baseflow,sigma_tab] = SF_FreeMovement_Spectrum('modefollow',baseflow,sigma_tab,RealShift,ImagShift,Stiffness_table,mass,nev);
    filename={[modename{1} '_data']}; %For the saved data (it's a cell)
    SF_Save_Data('data',General_data_dir,savedata_dir,Re,m_star,filename,Stiffness_table,U_star,sigma_tab);
    close all
end

%% Data Treatement: EigenValue 


%CHOOSE Data to plot: (Be sure that data exists)
%NOT FULLY OPERATIONAL YET
General_data_dir_folder=General_data_dir; % e.g.: './FOLDER_TOTO/'
domian_plot=domain_identity;        %e.g.:{'totodir1','totodir2'}
mesh_plot=mesh_identity;
%folder_plot={[General_data_dir_folder 'Post-Treatement/' domian_plot{1} mesh_plot{1} ]}; % isto funciona se forem 2 meshs pe?? 
folder_plot={[General_data_dir_folder  domian_plot{1} mesh_plot{1} ]}; % isto funciona se forem 2 meshs pe?? 

Re_plot=[19 19.8 19.9 20 23] ; % Re for previous calculation; for an array: [Re1 Re2]
m_star_plot=[20]; % m_star for previous calculation

%The different data treatment options:
%Mode:Fluid, Structure or Both
%Axis:sigma_VS_Ustar, F_LSA_VS_Ustar, f_star_LSA_VS_Ustar, sigma_r_VS_Ustar_LSA,
%NavroseMittal2016LockInRe60M20, NavroseMittal2016LockInRe60M5,
%NavroseMittal2016LockInRe40M10

%How to do:'Mode:-', 'Axis:-'

SF_Data_Treatement('Mode:Structure','Axis:sigma_VS_Ustar',General_data_dir,folder_plot,Re_plot,m_star_plot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If multiple data is ploted, RUN these lines for saving:
filename={'IMAGE_TOTO5'};%Name of the figure if it's a comparing data figure:
SF_Save_Data('graphic',General_data_dir,folder_plot,Re_plot,m_star_plot,filename,0,0,0); %Last 3 not used in 'graphic'

%% EigenMode Plot a faire
%CHOOSE the eigenmode to plot
%
%
%
%Computation of the demanding eigenmodes
%
%
%
%Eigenmode plot
%
%
%
%Eigenmode save data and figures
%chosen_mode=30;
%plotFF(mode_tab(chosen_mode),'vort1')





%% PUBELLE apartir d'ici:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Real and imaginary sigma values in function of Ustar
%SF_FreeMovement_Display(U_star,sigma_tab,'sigma_VS_Ustar',modename);
%filename=[ numbermode 'mode' modename];
%SF_Save_Data('grafic',savedata_dir,Re,m_star,filename,Stiffness_table,U_star,sigma_tab);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Real sigma and F_LSA values in function of Ustar(iqual to Navrose p572)
%SF_FreeMovement_Display(U_star,sigma_tab,'F_LSA_VS_Ustar',modename);
%filename=[ numbermode 'mode' modename '_Navrosep572'];
%Save_Data(Re,m_star,Stiffness_table,U_star,filename,sigma_tab,savedata_dir_version,'grafic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Navrose page578 (eg:Re40) same units:
%SF_FreeMovement_Display(U_star,sigma_tab,'f_LSA_VS_Ustar',modename);
%filename=[ numbermode 'mode' modename '_Navrosep578'];
%Save_Data(Re,m_star,Stiffness_table,U_star,filename,sigma_tab,savedata_dir_version,'grafic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%After have done the two modes, we can continue to the next code part:
%CHOOSE the same modenames as saved in the folder:
% modename1='STRU'; % Options normally used: STRU, FEM2
% modename2='FLUI';% Options normally used: FLUI, FEM1
% modename=[modename1 ; modename2];
% sigma_tab_both=[];
%  
% 
% sigma_dir=['./Final_results_' savedata_dir_version '/Re' num2str(Re) '/' 'mstar' num2str(m_star) '/'];
% numbermode='02'; load([sigma_dir numbermode 'mode' modename1 '_spectrum.mat']);
% sigma_tab_both(1,:)=sigma_tab;
% numbermode='03'; load([sigma_dir numbermode 'mode' modename2 '_spectrum.mat']);
% sigma_tab_both(2,:)=sigma_tab;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%One image iqual to Navrose p.572:
%SF_FreeMovement_Display(U_star,sigma_tab_both,'F_LSA_VS_Ustar',modename);
%filename='04mode_both';
%Save_Data(Re,m_star,Stiffness_table,U_star,filename,sigma_tab,savedata_dir_version,'grafic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image de Navrose page574
SF_FreeMovement_Display(U_star,sigma_tab_both,'Ustar_LSA_VS_sigma_i',modename);
filename='05_Navrose_snake';
%Save_Data(Re,m_star,Stiffness_table,U_star,filename,sigma_tab,savedata_dir_version,'grafic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Navrose page578 (eg:Re40) same units:
SF_FreeMovement_Display(U_star,sigma_tab_both,'f_LSA_VS_Ustar',modename);
filename='06_Navrose_p578';
%Save_Data(Re,m_star,Stiffness_table,U_star,filename,sigma_tab,savedata_dir_version,'grafic');

%End

%%%% Data Treatement of the values we want IMAGEM RE40 MASS VARIATION
%First, ensure that the are already calculated
%savedata_dir_version='v13';
Re_disp=[40];
modename='STRU'; % Options normally used: STRU, FEM2...
numbermode='02';
m_star_disp=[20 10 5 0.7];
i=1;
sigma_tab_both=[];

for Re=Re_disp
    for m_star=m_star_disp
    sigma_dir=['./Final_results_' savedata_dir_version '/Re' num2str(Re) '/' 'mstar' num2str(m_star) '/'];
    load([sigma_dir numbermode 'mode' modename '_spectrum.mat']);
    sigma_tab_both(i,:)=sigma_tab;
    i=i+1;
    end
end
SF_FreeMovement_Display(U_star,sigma_tab_both,'F_LSA_VS_Ustar',modename);
title('Oscillation rate (mass variation for Re40)')
currentFigure = gcf;
title(currentFigure.Children(end), 'Amplification rate (mass variation for Re40)');
%
legend('mstar=20','mstar=10','mstar=5','mstar=0.7')

saveas(gcf,['./Final_results_v13/massvariationRe40.fig']);
saveas(gcf,['./Final_results_v13/massvariationRe40.png']);

%%%% Data Treatement of the values we want IMAGEM RE VARIATION
%First, ensure that the are already calculated
%savedata_dir_version='v13';
Re_disp=[15 19 21];
modename='STRU'; % Options normally used: STRU, FEM2...
numbermode='02';
m_star_disp=[100];
i=1;
sigma_tab_both=[];

for Re=Re_disp
    for m_star=m_star_disp
    sigma_dir=['./Final_results_' savedata_dir_version '/Re' num2str(Re) '/' 'mstar' num2str(m_star) '/'];
    load([sigma_dir numbermode 'mode' modename '_spectrum.mat']);
    sigma_tab_both(i,:)=sigma_tab;
    i=i+1;
    end
end
SF_FreeMovement_Display(U_star,sigma_tab_both,'F_LSA_VS_Ustar',modename);
title('Oscillation rate (Re variation for mstar=100)')
currentFigure = gcf;
title(currentFigure.Children(end), 'Amplification rate (Re variation for mstar=100)');
%
legend('Re=15','Re=19','Re=21')

saveas(gcf,['./Final_results_v13/Revariation_mstar100.fig']);
saveas(gcf,['./Final_results_v13/Revariation_mstar100.png']);


%% Reynolds Analysis: Search for the critic Re curve
%U_star=[3:0.2:6 6:0.05:9.5 9.5:0.2:11];
U_star=[3:0.2:7 7:0.05:9.5 9.5:0.2:11];
%U_star=[3:0.2:11];
m_star_tab=[100:-25:50 50:-10:20 15:-5:5];
Re_init=47;
Re_before=47;
shift_init_struc=0+2*1i;
shift_init_fluid=0.02+0.75*1i;
stable=false;
table_plot=[];

figure(50); hold on;
for m_star=m_star_tab
    mass=m_star*pi/4;
    Stiffness_table=pi^3*m_star./((U_star).^2);
    
    Re_tab = Re_init: -0.5 : 30;
    for Re=Re_tab
        
        %FAIRE LE BASEFLOW
        
        [stable,shift_init_struc,shift_init_fluid]= SF_Stab_Curve_function(baseflow,shift_init_struc,shift_init_fluid,Stiffness_table,mass);
        
        if(stable)%==true
            table_plot=[table_plot m_star;Re];
            scatter(50,Re,m_star);
            Re_init=Re_before;
            break;
        end
        Re_before=Re;
    end
end

figure; hold on;
plot(table_plot(1),table_plot(2));


%% Trash
if(1==0)
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

%% starting point ; autres choses 


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
