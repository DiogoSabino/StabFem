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
plotFF(baseflow,'mesh');%pause(0.1);
baseflow=SF_BaseFlow(baseflow,'Re',60);
% 
disp(' ');disp('ADAPTING MESH FOR RE=60 ACORDING TO EIGENVALUE ');disp(' ');
%Mesh adaptation to a fluid mode
[ev,em] = SF_Stability(baseflow,'shift',0.04+0.74i,'nev',1,'type','D');


[baseflow,em]=SF_Adapt(baseflow,em,'Hmax',10,'InterpError',0.005);
plotFF(baseflow,'mesh');%pause(0.1);


% path of the saved data for the harmonic case cylinder (repeated in macros.edp)
ffdataharmonicdir=[ffdatadir 'DATA_SF_CYLINDER/'];


%% Erase previous data if mesh have change (do it manually!!)
% Delete the data concerning the solution of the problem
% (not the baseflow and the mesh; the others eg DATA_SF_CYLINDER)
% Careful in the use of the next lines
while(false) % change to 'true' if you want delete all this data
    system(['rm -R ' ffdataharmonicdir]);
    
    break % Compulsory to exit the while loop
end


%% Case of Harmonic Forcing Cylinder
% NB:verbosity= 10;if simulation seems to stuck do Ctrl+C and start again;
% previous data will not be lost and script will continue where it was interrupted
%clear all
%close all

%Data
if(exist(ffdataharmonicdir)~=7&&exist(ffdataharmonicdir)~=5)
    mysystem(['mkdir ' ffdataharmonicdir]);
end %It's compulsory: Creation of the directory

%Re_tab=20 %for testing just one value Re=20
Re_tab=[21]; %for testing multiple values % first run
%Re_tab=[19.5 19.6 19.7 19.8 19.9 20]; %For Re1 first try
%i didnt use it%Re_tab=[19.8 19.9 19.91 19.92 19.95 20]; %For Re1 SECOND try
%Re_tab=[30 30.1 30.2 30.3 30.4 30.5 32 ]% search for the 2nd crit Re2
%Re_tab=[46 46.5 46.6 46.7 46.8 46.9 47]% search for the 3nd crit Re3
%i didnt use it%Re_tab=[19.91 30.3] curves of the 3 Rec

Omega_values=[0.1:0.02:1.2]; % First Value; ...intermediate values ... ; Last value

verbosity=10; %To see all

for Re=Re_tab
   baseflow=SF_BaseFlow(baseflow,'Re',Re); 
   SF_HarmonicForcing(baseflow,Omega_values);
   disp('___________________________________________________')%To separete iters
end

%% Data treatement
%select the desired omegas and re

figure
Re_tab=[21];
for Re=Re_tab
    %extract data from importFFdata...
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    
    SF_HarmonicForcing_Display(dataRe.OMEGAtab,dataRe.Lift,Re_tab); 
end

%% Iterative method to find Rec1
%St_tab=Omega_values/2/pi;
%Omega_values=St_tab*2*pi;

Omega_values=0.65:0.005:0.67;
Re=19.9; %Re initial

baseflow=SF_BaseFlow(baseflow,'Re',Re);
SF_HarmonicForcing(baseflow,Omega_values);
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe=importFFdata(all_data_stored_file);
Zr_min=min(real(dataRe.Lift));

increment=0.1;
Re_last=25;% just to enter the loop
Zr_min_last=Zr_min;
Rec1_convergence_tab=[Re];


while(abs(Re_last-Re)>0.001)
    disp('___________________________________________________')%To separete iters
    Zr_min_last=Zr_min;
    Re_last=Re;
    
    if(Zr_min_last>0)
        Re=Re_last+increment;
        disp('TOTO1')
    elseif(Zr_min_last<0)
        Re=Re_last-increment;
        disp('TOTO2')
    elseif(Zr_min_last==0)
        disp('TOTO3')
        break;
    end
    baseflow=SF_BaseFlow(baseflow,'Re',Re);
    SF_HarmonicForcing(baseflow,Omega_values);
    
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    Zr_min=min(real(dataRe.Lift));
    
    if((Zr_min_last*Zr_min)<0)
        disp('TOTO4')
        increment=0.5*increment;
    end
    Rec1_convergence_tab=[Rec1_convergence_tab Re];
    disp(['Re=' num2str(Re)]);
    disp(['Re_last=' num2str(Re_last)]);
end

disp('Loop Terminated');
disp(['Rec1=' num2str(Re)]);
%19.9461

%% Iterative method to find Rec2
%St_tab=Omega_values/2/pi;
%Omega_values=St_tab*2*pi;

Omega_values=0.57:0.005:0.61;
Re=30.2; %Re initial

baseflow=SF_BaseFlow(baseflow,'Re',Re);
SF_HarmonicForcing(baseflow,Omega_values);
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe=importFFdata(all_data_stored_file);
Z_min=min(imag(dataRe.Lift));

increment=0.1;
Re_last=32;% just to enter the loop
Z_min_last=Z_min;
Rec2_convergence_tab=[Re];


while(abs(Re_last-Re)>0.001)
    disp('___________________________________________________')%To separete iters
    Z_min_last=Z_min;
    Re_last=Re;

    
    if(Z_min_last>0)
        Re=Re_last+increment;
    elseif(Z_min_last<0)
        Re=Re_last-increment;
    elseif(Z_min_last==0)
        break;
    end
    baseflow=SF_BaseFlow(baseflow,'Re',Re);
    SF_HarmonicForcing(baseflow,Omega_values);
    
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    Z_min=min(imag(dataRe.Lift));
    
    
    if((Z_min_last*Z_min)<0)
        increment=0.5*increment;
    end
    
    Rec2_convergence_tab=[Rec2_convergence_tab Re];
    disp('End of iteration')
    disp(['Re_last=' num2str(Re_last)]);
    disp(['Re=' num2str(Re)]);
    disp(['min Zi=' num2str(Z_min)]);
    disp(['increment for next iter=' num2str(increment)]);

end

disp('Loop Terminated');
disp(['Rec2=' num2str(Re)]);

%% Iterative method to find Rec3: not the best one, but it works...
%St_tab=Omega_values/2/pi;
%Omega_values=St_tab*2*pi;

%The refinement of omega is extremly important in this iteration method...
%The best would be a automatic refinement of omega...
Omega_values=0.73:0.0005:0.735;
Re=46.6; %Re initial

baseflow=SF_BaseFlow(baseflow,'Re',Re);
SF_HarmonicForcing(baseflow,Omega_values);
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe=importFFdata(all_data_stored_file);

Yr_close=min(abs(real(1./dataRe.Lift)));
index_close=find(abs(Yr_close-abs(real(1./dataRe.Lift)))<10^(-6) );
Yi_close=imag(1/dataRe.Lift(index_close));

norm=sqrt( real(1/dataRe.Lift(index_close)).^2+ imag(1/dataRe.Lift(index_close))^2 ) ;
norm_tab=[norm];
increment=0.1;
Re_last=46;% just to enter the loop
Yi_close_last=Yi_close;
Rec3_convergence_tab=[Re];


while(abs(Re_last-Re)>0.001)
    disp('___________________________________________________')%To separete iters
    Yi_close_last=Yi_close;
    Re_last=Re;
    
    if(Yi_close>0)
        Re=Re_last+increment;
        disp('TOTO1')
    elseif(Yi_close<0)
        Re=Re_last-increment;
        disp('TOTO2')
    elseif(Yi_close==0)
        disp('TOTO3')
        break;
    end
    baseflow=SF_BaseFlow(baseflow,'Re',Re);
    SF_HarmonicForcing(baseflow,Omega_values);
    
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    
    Yr_close=min(abs(real(1./dataRe.Lift)));
    index_close=find(abs(Yr_close-abs(real(1./dataRe.Lift)))<10^(-6) );
    Yi_close=imag(1/dataRe.Lift(index_close));
    
    norm=sqrt( real(1/dataRe.Lift(index_close)).^2+ imag(1/dataRe.Lift(index_close))^2 ) ;
    
    if((Yi_close_last*Yi_close)<0)
        disp('TOTO4')
        increment=0.5*increment;
    end
    Rec3_convergence_tab=[Rec3_convergence_tab Re];
    norm_tab=[norm_tab norm];
    disp(['Re=' num2str(Re)]);
    disp(['Re_last=' num2str(Re_last)]);
end

disp('Loop Terminated');
disp(['Rec3=' num2str(Re)]);

%% Impedance Predictions
%See IMPEDANCE_based_predictions.m

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







