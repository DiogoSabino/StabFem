%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is part of StabFem Project.
%   (set of programs in Freefem to perform global stability calculations  
%   in hydrodynamic stability)
%
%	File: SCRIPT_VIV_HARMONICFORCED.m
%   Contributours: David Fabre, Diogo Sabino ...
%   Last Modification: Diogo Sabino, 18 April 2018
%   
%	This script models a flow aroung a rigid circular cilinder mounted
%	in a spring-damped system to simulate VIV
%
%   To give a good use of this script, run each section at a time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh creation and convergence using adapt mesh
%clear all
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

%% Erase previous data if mesh have change (do it manually!!)
% Delete the data concerning the solution of the problem
% (not the baseflow and the mesh; the others eg DATA_SF_CYLINDER)
% NB: If one wants to delete just some of the values of .txt file, go do it manualy

% path of the saved data for the harmonic case cylinder (repeated in macros.edp)
ffdataharmonicdir=[ffdatadir 'DATA_SF_CYLINDER/'];

% Careful in the use of the next lines
while(false) % change to 'true' if you want delete all this data
    system(['rm -R ' ffdataharmonicdir]);
    
    break % Compulsory to exit the while loop
end


%% Case of Harmonic Forcing Cylinder
% NB:verbosity= 10;if simulation seems to stuck do Ctrl+C and start again;
% previous data will not be lost and script will continue where it was interrupted
%clear all
%disp('Press any key to continue to Harmonic Forcing Case');
%pause; 
%close all

%Data
if(exist(ffdataharmonicdir)~=7&&exist(ffdataharmonicdir)~=5) %je n'est pas compris tres bien cette commande; a voir ensemble apres
    mysystem(['mkdir ' ffdataharmonicdir]);
end %It's compulsory: Creation of the directory

%Re_tab=20 %for testing just one value Re=20
Re_tab=[15 20:10:60]; %for testing multiple values % first run
%Re_tab=[19.5 19.6 19.7 19.8 19.9 20]; %For Re1 first try
%i didnt use it%Re_tab=[19.8 19.9 19.91 19.92 19.95 20]; %For Re1 SECOND try
%Re_tab=[30 30.1 30.2 30.3 30.4 30.5 32 ]% search for the 2nd crit Re2
%Re_tab=[46 46.5 46.6 46.7 46.8 46.9 47]% search for the 3nd crit Re3
%i didnt use it%Re_tab=[19.91 30.3] curves of the 3 Rec

%NB: If one adds a Re and wants to have the same data of the other Re already
%       calculated... well, this isn't code yet, the best solution is to 
%       delete all data and computed it again for all the Re needed


%Don't give silly values, sometimes matlab inexplicably crash with them
% First Value; Step; Last value
Omega_values=0.3:0.02:1.2; % First Value;Step;Last value|eg 0.3:0.02:1.2

Omega_fine=0; %First run
%Omega_fine=0.62:0.005:0.68;   %Refinement for Rec1
%Omega_fine=0.59:0.0005:0.61;  %Refinement for Rec2
%Omega_fine=0.65:0.005:0.75;  %Refinement for Rec3
%
if isempty(Omega_fine)==1 %Create a array of omega_values to test
    disp('Minor error: omega raffinement not achieved');
elseif isempty(Omega_values)==1
    disp('Major error: Omega values missed!');
    Omega_values=-1;
else
    for i = 1:numel(Omega_values);
        if Omega_values(i) >= Omega_fine(1)
            n_upper=find(Omega_values>Omega_fine(end));
            Omega_values=[Omega_values(1:i-1) Omega_fine Omega_values(n_upper:end)];
            break;
        end
    end
end
%Omega_values will be between 
%min(Omega_values,Omega_fine) and max(Omega_values,Omega_fine)

verbosity=10; %To see all


figure
for Re=Re_tab
   baseflow=SF_BaseFlow(baseflow,'Re',Re); 
   [Omegatab,Ltab,Mtab]=SF_HarmonicForcing(baseflow,Omega_values);
   SF_HarmonicForcing_Display(Omegatab,Ltab,Re_tab);
   disp('___________________________________________________')%To separete iters
end


%% Spectial treatement for a beautiful plot of Rec3: do it manually
if(1==0)
    %First erase the purple line manually in graph 5 for Re=46.7
    %The execute this
    txt_data = importdata([ ffdataharmonicdir,'HARMONIC_2D_LateralOscillations_Re46.7.txt']);
	Omegatab= txt_data(:,1);
	Ltab = txt_data(:,3)+1i*txt_data(:,4);
    U=1; D=1; %Cylinder non-dimentionnal definitions
    Strouhal_number=Omegatab/(2*pi)*U/D;
    subplot(2,3,5)
    plot(real(Ltab(1:35)),imag(Ltab(1:35)),'m-o','MarkerSize',2);
    plot(real(Ltab(36:end)),imag(Ltab(36:end)),'m-o','MarkerSize',2);
    %Finally change manually the collor of the line to correspond to the
    %legend. Right-click mouse-> edit...etc ;)
end

