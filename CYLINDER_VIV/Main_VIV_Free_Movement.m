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
global  verbosity ffdatadir

%CHOOSE the domain parameters:
domain_parameters=[-50 50 50];

[baseflow]= SF_MeshGeneration(domain_parameters);

%% Save Data

%CHOOSE folder for saving data:
General_data_dir='./Final_Results_v20/'; %a array of char

domain_identity={[ num2str(domain_parameters(1)) '_' num2str(domain_parameters(2)) '_' num2str(domain_parameters(3)) '/']};
%CHOOSE the name of the folder mesh: (It's good to choose the name because we can adapt mesh during)
mesh_identity={'Adapt_mode_Hmax1_InterError_0.02/'};
savedata_dir={[ General_data_dir domain_identity{1} mesh_identity{1}]};
%savedata_dir=[General_data_folder domain_identity mesh_identity];

%% Computation: Spectrum
%Parameters' Definition
%CHOOSE the Re to test:
Re_tab=[40]; verbosity=10;
for Re=Re_tab
    baseflow = SF_BaseFlow(baseflow,'Re',Re);
    
    %CHOOSE the m_star to test:
    m_star_tab=[70];
    
    for m_star=m_star_tab
        mass=m_star*pi/4;
        %CHOOSE U_star
        
        %U_star=[1:0.05:4]; %m=0.05
        %U_star=[2:0.1:7]; %m=0.2,0.1
        %U_star=[3:0.1:7]; %m=1...0.8
        U_star=[3 3.05 3.1 3.2:0.2:6.5 6.5:0.05:11];
        %U_star=[11.1:0.2:20];%
        %U_star=[20:0.2:28];%
        
        Stiffness_table=pi^3*m_star./((U_star).^2);
        
        % Spectrum search: See spectrum if wanted, for discover the shift
        if(1==0)
            sigma_tab=[];
            %Re and mass are already defined above. Define numbers or tables for:
            m=1000; % mass
            STIFFNESS_to_search=1000; %Use whatever we want, e.g.: Stiffness_table(1)
            RealShift=[0.05 0 -0.03];
            ImagShift=[0.75];
            
            nev=10; %Normally here we don't use just one
            [baseflow,sigma_tab] = SF_FreeMovement_Spectrum('search',baseflow,sigma_tab,RealShift,ImagShift,STIFFNESS_to_search,0,nev);
            filename={'01spectrum_search'};
            %SF_Save_Data('spectrum',General_data_dir,savedata_dir,Re,m_star,filename,Stiffness_table,U_star,sigma_tab);
        end
        
        % Mode Follow: Follow a mode along the Stiffness_table/U_star       
        %CHOOSE the one for save data w/ a good name:
        modename={'02modeSTRUCTURE'};
        %modename={'03modeFLUID'};
        
        [RealShift, ImagShift]=SF_Shift_selection(modename,Re,m_star);
        sigma_tab=[]; %RealShift=0.005; ImagShift=0.64;
        nev=1; %Normally if's just one, but if shift is wrong, it helps put more
        
        [baseflow,sigma_tab] = SF_FreeMovement_Spectrum('modefollow',baseflow,sigma_tab,RealShift,ImagShift,Stiffness_table,mass,nev);
        filename={[modename{1} '_data']}; %For the saved data (it's a cell)
        SF_Save_Data('data',General_data_dir,savedata_dir,Re,m_star,filename,Stiffness_table,U_star,sigma_tab);
        close all
    end
end

%% EigenValue: Data Treatement 

%CHOOSE Data to plot: (Be sure that data exists)
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1','totodir2'} %%FALTA POR O DOMAIN NO PROXIMO LOOP
mesh_plot={'Adapt_mode_Hmax10_InterError_0.02/','Adapt_mode_Hmax1_InterError_0.02/'};
folder_plot={};
for element=1:size(mesh_plot,2)
folder_plot{end+1}=[General_data_dir_folder  domain_plot{1} mesh_plot{element} ];
end

Re_plot=[45] ; % Re for previous calculation; for an array: [Re1 Re2]
m_star_plot=[70]; % m_star for previous calculation

%The different data treatment options: %How to do:'Mode:-', 'Axis:-'
%Mode:Fluid, Structure or Both
%Axis:sigma_VS_Ustar, F_LSA_VS_Ustar, f_star_LSA_VS_Ustar, sigma_r_VS_Ustar_LSA,
%NavroseMittal2016LockInRe60M20, NavroseMittal2016LockInRe60M5,
%NavroseMittal2016LockInRe40M10, spectrum, sigma_rCOMP,sigma_rCOMPRe33m50
%For testing:
%SF_Data_Treatement('Mode:Both','Axis:NavroseMittal2016LockInRe60M20',folder_plot,Re_plot,m_star_plot);

%IF one data only is plotted
filename=SF_Data_Treatement('Mode:Structure','Axis:spectrum',folder_plot,Re_plot,m_star_plot);


SF_Save_Data('graphic',General_data_dir_folder,folder_plot,Re_plot,m_star_plot,filename,0,0,0); %Last 3 not used in 'graphic'
%ELSE
SF_Data_Treatement('Mode:Both','Axis:sigma_rCOMP',folder_plot,Re_plot,m_star_plot);
filename={'Re60_m20_sacar_parte_real'};%CHOOSE name of the figure
SF_Save_Data('graphic',General_data_dir_folder,folder_plot,Re_plot,m_star_plot,filename,0,0,0); %Last 3 not used in 'graphic'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EigenMode: Plot
%CHOOSE the eigenmode to plot
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1','totodir2'}
mesh_plot={'Adapt_mode_Hmax10_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]}; % isto funciona se forem 2 meshs ??

Re_plot=[60] ; % Re for previous calculation; for an array: [Re1 Re2]
m_star_plot=[20]; % m_star for previous calculation

%Usage:
%%% availability, Compute
%%%Mode:Fluid, Structure or Both

%SEE U_star available
SF_Mode_Display('availability','Mode:Fluid',0,folder_plot,Re_plot,m_star_plot,0);

%COMPUTE of the demanding eigenmodes
U_star_plot=[8]; %put just one for now...
baseflow = SF_BaseFlow(baseflow,'Re',Re_plot);
[em]=SF_Mode_Display('Compute','Mode:Fluid',baseflow,folder_plot,Re_plot,m_star_plot,U_star_plot);

%Eigenmode plot reffinement
%...to do...
plotFF(em,'ux1.re')
plotFF(em,'vort1')
plotFF(em,'vort1.re')
%exportFF_tecplot(em,'./Modes_Tecplot/Re60m20U8Fluid.plt')
%Baseflow
%plotFF(baseflow,'vort')
%% Grafic m* vs fequency: not working yet
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1','totodir2'}
mesh_plot={'Adapt_mode_Hmax10_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]}; % isto funciona se forem 2 meshs ??
Re_plot=[40];
m_star_plot=[10 20 50 75]; %for m_star we want
%USE 'Mode:Both','Axis:Zone_m_star_vs_Ustar'
SF_Data_Treatement('Mode:Both','Axis:Zone_m_star_vs_Ustar',folder_plot,Re_plot,m_star_plot);

%% Grafic U* Re
%Put just one path!!
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1','totodir2'}
mesh_plot={'Adapt_mode_Hmax10_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]}; % isto funciona se forem 2 meshs ??
Re_plot=[20 23 25 27 30 33 35 40 42 43 45];
m_star_plot=50;
%m_star_plot=4.73;

Point_TAB_ReUstar_plane=SF_Data_Tretement_Post('Case:Ustar_Re_plane',folder_plot,Re_plot,m_star_plot);
%save( './Latex_data/Free/Koucomp/m473_ReUstar_plane.mat','Point_TAB_ReUstar_plane');
%save( './Latex_data/Free/Koucomp/m50_ReUstar_plane.mat','Point_TAB_ReUstar_plane');


%% Impedance values
%CHOOSE Data to plot: (Be sure that data exists)
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1','totodir2'}
mesh_plot={'Adapt_mode_Hmax10_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]}; % isto funciona se forem 2 meshs ??

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rec1
Re_plot=[19.95] ; 
m_star_plot=[0.05 0.1 0.2 0.3 0.5 1 2 3 4.73 5 10 20];

Point_TAB_Rec1=SF_Data_Tretement_Post('Case:Impedance_Points',folder_plot,Re_plot,m_star_plot);
%save( './Impedance_Treatement/FreeCaseRec1.mat','Point_TAB_Rec1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% David curve
%Re=21
Re_plot=[21] ;
m_star_plot=[100];

extracting=[folder_plot{1} 'Re' num2str(Re) '/mstar' num2str(m_star_plot) '/02modeSTRUCTURE_data.mat'];
DataFree_curveUstar_Re21_m100=load(extracting,'sigma_tab','U_star');
%save( './Impedance_Treatement/FreeCaseRe21mstar100.mat','DataFree_curveUstar_Re21_m100');
%Re=40
Re_plot=[40] ;
m_star_plot=[100];

extracting=[folder_plot{1} 'Re' num2str(Re) '/mstar' num2str(m_star_plot) '/02modeSTRUCTURE_data.mat'];
DataFree_curveUstar_Re40_m100=load(extracting,'sigma_tab','U_star');
%save( './Impedance_Treatement/FreeCaseRe40mstar100.mat','DataFree_curveUstar_Re40_m100');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparacao com Kou
%... a fazer
    
%% Put good legends
if(1==0)
    
    legend('Mode Fluid: Re60,mstar20','Mode Structure: Re60,mstar20', ...
        'Navrose Mode Structure: Re60,mstar20','Navrose Mode Fluid: Re60,mstar20')
    legend('Mode Fluid: Re60,mstar5','Mode Structure: Re60,mstar5', ...
        'Navrose Mode Structure: Re60,mstar5','Navrose Mode Fluid: Re60,mstar5')
    legend('Mode Structure: Re60,mstar20','Navrose','Zhang')

      legend('Diogo: Mode Fluid','Diogo: Mode Structure','Navrose: Mode Structure',...
          'Navrose: Mode Fluid','Zhang: Mode Structure','Zhang: Mode Fluid')
          legend('Diogo: Mode Structure','Zhang: Mode Structure')
end

%% Critic Reynolds Analysis: Search for the critic Re curve

%To do FindThreshold_VIV.epd
%a faire
%[baseflowC,emC]=SF_FindThreshold_VIV(baseflow,em);

%% Non-linear: Harmonic Balance
%% first try

%Re near the threshold

Re_start=20.5;
baseflow = SF_BaseFlow(baseflow,'Re',Re_start);
m_star=300;
U_star=9.6;

mass=m_star*pi/4;
STIFFNESS=pi^3*m_star./((U_star).^2);

shift=0+0.655i;
[ev,em] = SF_Stability(baseflow,'shift',shift,'nev',1,'type','D','STIFFNESS',STIFFNESS,'MASS',mass,'DAMPING',0,'Frame','R');

YGuess=0.01;
[meanflow,mode] = SF_SelfConsistentDirect(baseflow,em,'Yguess',YGuess,'STIFFNESS', STIFFNESS,'MASS',mass,'DAMPING',0);

%% second try

Re_start=45;
baseflow = SF_BaseFlow(baseflow,'Re',Re_start);
m_star=70;
U_star=6.4;

mass=m_star*pi/4;
STIFFNESS=pi^3*m_star./((U_star).^2);

shift=0+1i;
[ev,em] = SF_Stability(baseflow,'shift',shift,'nev',1,'type','D','STIFFNESS',STIFFNESS,'MASS',mass,'DAMPING',0,'Frame','R');

YGuess=0.2;
Amplitude_HB=[];
[meanflow,mode] = SF_SelfConsistentDirect(baseflow,em,'Yguess',YGuess,'STIFFNESS', STIFFNESS,'MASS',mass,'DAMPING',0);

%%
%data=importFFdata(baseflow.mesh,[ffdatadir 'SelfConsistentMode.ff2m']);
%Amplitude_HB=[Amplitude_HB data.Y];
Amplitude_HB=[Amplitude_HB mode.Y];



Re_tab=[46 48 55 60 65 70 75 80 85 90 95 100];
%Re_tab=[65 70 75 80 85 90 95 100];

for Re=Re_tab
    [meanflow,mode]=SF_SelfConsistentDirect(meanflow,mode,'Re',Re,'STIFFNESS', STIFFNESS,'MASS',mass,'DAMPING',0);
    Amplitude_HB=[Amplitude_HB mode.Y];  
end
figure
plot([45 Re_tab],real(Amplitude_HB))

system(['cp ' ffdatadir 'SelfConsistentMode.ff2m ' ffdatadir 'HB_O1/HB_ModeO1_mstar70_Re100_Ustar6p4.ff2m'])
system(['cp ' ffdatadir 'MeanFlow.ff2m ' ffdatadir 'HB_O1/MeanFlow_mstar70_Re100_Ustar6p4.ff2m'])
%save('./WORK/HB_O1/amplitude_HBO1_untilRe100.mat','Amplitude_HB', 'Re_tab')


%% follow at Re=100

system(['cp ' ffdatadir 'HB_O1/HB_ModeO1_mstar70_Re100_Ustar6p4.ff2m ' ffdatadir 'SelfConsistentMode.ff2m ' ])
system(['cp ' ffdatadir 'HB_O1/MeanFlow_mstar70_Re100_Ustar6p4.ff2m ' ffdatadir 'MeanFlow.ff2m '])

Re=100;
m_star=70;
U_star_plot=[6.1 6 5.8 5.6 5.4 5.2];
Amplitude_varying_Ustar=Amplitude_HB(end);

for U_star=U_star_plot
    STIFFNESS=pi^3*m_star./((U_star).^2);
    [meanflow,mode]=SF_SelfConsistentDirect(meanflow,em,'Re',Re,'STIFFNESS', STIFFNESS,'MASS',mass,'DAMPING',0);
    Amplitude_varying_Ustar=[Amplitude_varying_Ustar mode.Y];
end
figure
plot([6.2 U_star_plot],real(Amplitude_varying_Ustar))

%save('./WORK/HB_O1/amplitude_HBO1_varyingUstar.mat','Amplitude_varying_Ustar','U_star_plot' )

    