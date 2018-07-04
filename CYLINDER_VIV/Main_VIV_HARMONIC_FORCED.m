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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh: creation and convergence using adapt mesh
%clear all
global ffdataharmonicdir verbosity ffdatadir
%CHOOSE the domain parameters:
domain_parameters=[-50 50 50];

[baseflow]= SF_MeshGeneration(domain_parameters);

%% path of the saved data for the harmonic case cylinder (repeated in macros.edp)
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

%%% First Run
Re_tab=[46.6];
Omega_values=[0.5:0.005:1.2];
%%% Second Run refining where needed (first figures)
%Re_tab=[25 35];
%Omega_values=[0.46:0.005:0.74];
%Re_tab=[55];
%Omega_values=[0.6:0.005:0.86]; %Additional values for Re=55
%%% Third run (Curves for 21 and 100)
%Re_tab=[21 40];
%Omega_values=[0.38:0.005:0.82];
%%% For critic Re
%Re_tab=20;
%Omega_values=[0.1:0.005:1.2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbosity=10; %To see all
for Re=Re_tab
   baseflow=SF_BaseFlow(baseflow,'Re',Re); 
   SF_HarmonicForcing(baseflow,Omega_values);
   disp('___________________________________________________')%To separete iters
    %Calculating derivative:
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    dZr=diff(2*real(dataRe.Lift));
    dZi=diff(2*real(dataRe.Lift));
    save([ffdataharmonicdir 'Forced_Re' num2str(Re) '_diff_Lift_Coeff.mat'],'dZr','dZi');
end

%% Data treatement
%select the desired omegas and re

figure
Re_tab=[46.6 46.7 46.7688 46.7656 46.7641 46.7648 46.7652 ];
for Re=Re_tab
    %extract data from importFFdata...
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    SF_HarmonicForcing_Display(dataRe.OMEGAtab,dataRe.Lift,Re_tab); 
end

%% Iterative method to find Rec1

%erase previous values for Re=20 for assured convergence
Re=20; %Re initial
Omega_values=linspace(0.64,0.67,10);%omegas init
baseflow=SF_BaseFlow(baseflow,'Re',Re);

SF_HarmonicForcing(baseflow,Omega_values);
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe=importFFdata(all_data_stored_file);
[Zr_min,indexmin]=min(real(dataRe.Lift));
Omegamin=dataRe.OMEGAtab(indexmin);
Zimin=imag(dataRe.Lift(indexmin));

incrementRe=0.1;
incrementOMEGA=0.02;
Re_last=25;% just to enter the loop
Zr_min_last=Zr_min;
Rec1_convergence_tab=[Re];
incrementOMEGA_TAB=[0.02];
Omegamin_tab=[Omegamin];
Zimin_tab=[Zimin];
indexmin_tab=[indexmin];

while(abs(Re_last-Re)>0.0005 || incrementOMEGA>0.0005)
    disp('___________________________________________________')%To separete iters
    Zr_min_last=Zr_min;
    Re_last=Re;
    
    if(Zr_min_last>0)
        Re=Re_last+incrementRe;
    elseif(Zr_min_last<0)
        Re=Re_last-incrementRe;
    elseif(Zr_min_last==0)
        break;
    end
    
    baseflow=SF_BaseFlow(baseflow,'Re',Re);
    Omega_values=linspace(Omegamin-incrementOMEGA,Omegamin+incrementOMEGA,10);
    SF_HarmonicForcing(baseflow,Omega_values);
    
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    [Zr_min,indexmin]=min(real(dataRe.Lift));
    Omegamin=dataRe.OMEGAtab(indexmin);
    
    if((Zr_min_last*Zr_min)<0&& abs(Re_last-Re)>0.0005)
        incrementRe=0.5*incrementRe;
    end
    if((Zr_min_last*Zr_min)<0 && incrementOMEGA>0.0005)
        incrementOMEGA=0.5*incrementOMEGA;
    end
    
    Rec1_convergence_tab=[Rec1_convergence_tab Re];
    incrementOMEGA_TAB=[incrementOMEGA_TAB incrementOMEGA];
    Omegamin_tab=[Omegamin_tab Omegamin];
    Zimin_tab=[Zimin_tab Zimin];
    indexmin_tab=[indexmin_tab indexmin];
    disp(['RE TAB:' num2str(Rec1_convergence_tab)]);
end

disp('Loop Terminated');
disp(['Rec1=' num2str(Re)]);
%Re=19.9152
%omega=0.6533
%St=0.1040
%Zi=0.9504*2=1.9010

%% Iterative method to find Rec2

%erase previous values for Re=30.3 for assured convergence
Re=30.3; %Re initial
Omega_values=linspace(0.56,0.62,10);%omegas init
baseflow=SF_BaseFlow(baseflow,'Re',Re);

SF_HarmonicForcing(baseflow,Omega_values);
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe=importFFdata(all_data_stored_file);
[Zi_min,indexmin]=min(imag(dataRe.Lift));
Omegamin=dataRe.OMEGAtab(indexmin);
Zrmin=imag(dataRe.Lift(indexmin));

incrementRe=0.1;
incrementOMEGA=0.02;
Re_last=32;% just to enter the loop
Zi_min_last=Zi_min;
Rec2_convergence_tab=[Re];
incrementOMEGA_TAB=[0.02];
Omegamin_tab=[Omegamin];
Zrmin_tab=[Zrmin];
indexmin_tab=[indexmin];

while(abs(Re_last-Re)>0.0005 || incrementOMEGA>0.0005)
    disp('___________________________________________________')%To separete iters
    Zi_min_last=Zi_min;
    Re_last=Re;
    
    if(Zi_min_last>0)
        Re=Re_last+incrementRe;
    elseif(Zi_min_last<0)
        Re=Re_last-incrementRe;
    elseif(Zi_min_last==0)
        break;
    end
    
    baseflow=SF_BaseFlow(baseflow,'Re',Re);
    Omega_values=linspace(Omegamin-incrementOMEGA,Omegamin+incrementOMEGA,10);
    SF_HarmonicForcing(baseflow,Omega_values);
    
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    [Zi_min,indexmin]=min(imag(dataRe.Lift));
    Omegamin=dataRe.OMEGAtab(indexmin);
    
    if((Zi_min_last*Zi_min)<0&& abs(Re_last-Re)>0.0005)
        incrementRe=0.5*incrementRe;
    end
    if((Zi_min_last*Zi_min)<0 && incrementOMEGA>0.0005)
        incrementOMEGA=0.5*incrementOMEGA;
    end
    
    Rec2_convergence_tab=[Rec2_convergence_tab Re];
    incrementOMEGA_TAB=[incrementOMEGA_TAB incrementOMEGA];
    Omegamin_tab=[Omegamin_tab Omegamin];
    Zrmin_tab=[Zrmin_tab Zrmin];
    indexmin_tab=[indexmin_tab indexmin];
end

disp(['RE TAB:' num2str(Rec2_convergence_tab)]);
disp('Loop Terminated');
disp(['Rec2=' num2str(Re)]);
%Re=30.349
%omega=0.5972
%St=0.0950
%Zr=0

%St_tab=Omega_values/2/pi;
%Omega_values=St_tab*2*pi;

%% Iterative method to find Rec3
%St_tab=Omega_values/2/pi;
%Omega_values=St_tab*2*pi;

%erase previous values for Re=46.6 for assured convergence
Re=46.6; %Re init
Omega_values=linspace(0.72,0.75,10);%omegas init

baseflow=SF_BaseFlow(baseflow,'Re',Re);
SF_HarmonicForcing(baseflow,Omega_values);
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe=importFFdata(all_data_stored_file);

[Yr_close,index_close]=min(abs(real(1./dataRe.Lift)));
Yi_close=imag(1/dataRe.Lift(index_close));
Omegamin=dataRe.OMEGAtab(index_close);

norm=sqrt( real(1/dataRe.Lift(index_close)).^2+ imag(1/dataRe.Lift(index_close))^2 ) ;
norm_tab=[norm];
incrementRe=0.1;
incrementOMEGA=0.02;
Re_last=46;% just to enter the loop
Rec3_convergence_tab=[Re];
incrementOMEGA_TAB=[0.02];
Omegamin_tab=[Omegamin];
ALL_OMEGA_VALUES=[Omega_values];


while(abs(Re_last-Re)>0.0005 || incrementOMEGA>0.0005)
    disp('___________________________________________________')%To separete iters
    Yi_close_last=Yi_close;
    Re_last=Re;
    
    if(Yi_close>0)
        Re=Re_last+incrementRe;
    elseif(Yi_close<0)
        Re=Re_last-incrementRe;
    elseif(Yi_close==0)
        break;
    end
    
    baseflow=SF_BaseFlow(baseflow,'Re',Re);
    Omega_values=linspace(Omegamin-incrementOMEGA,Omegamin+incrementOMEGA,10);
    ALL_OMEGA_VALUES=[ALL_OMEGA_VALUES;Omega_values];
    SF_HarmonicForcing(baseflow,Omega_values);
    
    all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
    dataRe=importFFdata(all_data_stored_file);
    
    [Yr_close,index_close]=min(abs(real(1./dataRe.Lift)));
    Yi_close=imag(1/dataRe.Lift(index_close));
     Omegamin=dataRe.OMEGAtab(index_close);
    
    norm=sqrt( real(1/dataRe.Lift(index_close)).^2+ imag(1/dataRe.Lift(index_close))^2 ) ;
    
    if((Yi_close_last*Yi_close)<0&& abs(Re_last-Re)>0.0005)
        incrementRe=0.5*incrementRe;
    end
    if((Zi_min_last*Zi_min)<0 && incrementOMEGA>0.0005)
        incrementOMEGA=0.5*incrementOMEGA;
    end
    
    
    Rec3_convergence_tab=[Rec3_convergence_tab Re];
    norm_tab=[norm_tab norm];
    incrementOMEGA_TAB=[incrementOMEGA_TAB incrementOMEGA];
    Omegamin_tab=[Omegamin_tab Omegamin];
    
end

disp(['RE TAB:' num2str(Rec3_convergence_tab)]);
disp('Loop Terminated');
disp(['Rec3=' num2str(Re)]);
%Re=46.766 (see manually after covergence, due to lack of precision on omega)
%omega=0.7323
%St=0.1165
%Zr=0

%% TRASH: Spectial treatement for a beautiful plot of Rec3: do it manually
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
