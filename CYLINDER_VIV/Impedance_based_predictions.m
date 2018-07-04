%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Impedance method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Limit Rec1 %A REFAZER

global ffdataharmonicdir  ffdatadir

%run('../SOURCES_MATLAB/SF_Start.m');
%ffdataharmonicdir=[ffdatadir 'DATA_SF_CYLINDER/'];



Re=19.9461; %Rec1
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe_total=importFFdata(all_data_stored_file);


[Zr,index]=min(real(dataRe_total.Lift));
Zi=imag(dataRe_total.Lift(index));
Omega_f=dataRe_total.OMEGAtab(index);

%Variation of mstar:
mstar=[0:0.05:19.9 20];

Ustar_impedance=sqrt((4*(pi^2))./(Omega_f.*(Omega_f+2./(pi.*mstar).*Zi)));

%Plot predictions:
%figure
%plot(mstar,Ustar_impedance);
%figure
%loglog(mstar,Ustar_impedance);
%figure
%semilogx(mstar,Ustar_impedance);

% Load data from Free case:
FreeCase=load('./Impedance_Treatement/FreeCase.mat');
m_free=FreeCase.Point_TAB_Rec1(2,:);
U_free=FreeCase.Point_TAB_Rec1(3,:);

%Figure
figure;
plot(mstar,Ustar_impedance); hold on;
scatter(m_free,U_free,'r');

title('Predictions of Impedance Based Method');
xlabel('m^*'); ylabel('U^*');
legend('Impedance Based Predictions','Data from free Case')

%Figure with x-axis in log scale
figure;
semilogx(mstar,Ustar_impedance); hold on;
scatter(m_free,U_free,'r');

title('Predictions of Impedance Based Method (in log scale)');
xlabel('m^*'); ylabel('U^*');
legend('Impedance Based Predictions','Data from free Case')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OLD:Courbe que David demande...
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1'}
mesh_plot={'Adapt_mode_Hmax10_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Re=21
% For m=100
Re=21;
mstar=100;

%Load Data Forced-Case:
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe_total=importFFdata(all_data_stored_file);
%Data Forced case:
Zr=2*real(dataRe_total.Lift);
Zi=2*imag(dataRe_total.Lift);
Omega_f=dataRe_total.OMEGAtab;

Ustar_impedance=sqrt((4*(pi^2))./(Omega_f.*(Omega_f+2./(pi.*mstar).*Zi)));
lambda_r_impedance=-Zr./(mstar*pi);

%Load Data Free-Case:
FreeCase=load([folder_plot{1} 'Re' num2str(Re) '/mstar' num2str(mstar) '/02modeSTRUCTURE_data.mat']);
lambda_r_Free=real(FreeCase.sigma_tab);
U_free=FreeCase.U_star;

%Figures:
figure, hold on
plot(Ustar_impedance,lambda_r_impedance); %0.5 é por causa da differenca L e CL
plot(U_free,lambda_r_Free);
plot([8.5 10.5],[0 0],'k--')
title('Comparision Impedance-Based-Method and Free Results for Re=21');
xlabel('U^*'); ylabel('\lambda_r or Z_r');
legend('Impedance Based Predictions','Data from free Case')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Re=40
% For m=100
Re=40;
mstar=100;

%Load Data Forced-Case:
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe_total=importFFdata(all_data_stored_file);
%Data Forced case:
Zr=2*real(dataRe_total.Lift);
Zi=2*imag(dataRe_total.Lift);
Omega_f=dataRe_total.OMEGAtab;

Ustar_impedance=sqrt((4*(pi^2))./(Omega_f.*(Omega_f+2./(pi.*mstar).*Zi)));
lambda_r_impedance=-Zr./(mstar*pi);

% Load data from Free case:
FreeCase=load([folder_plot{1} 'Re' num2str(Re) '/mstar' num2str(mstar) '/02modeSTRUCTURE_data.mat']);
lambda_r_Free=real(FreeCase.sigma_tab);
U_free=FreeCase.U_star;

%Figures:
figure, hold on
plot(Ustar_impedance,lambda_r_impedance); %0.5 é por causa da differenca L e CL
plot(U_free,lambda_r_Free);
plot([5 20],[0 0],'k--')
title('Comparision Impedance-Based-Method and Free Results for Re=40');
xlabel('U^*'); ylabel('\lambda_r or Z_r');
legend('Impedance Based Predictions','Data from free Case')

%% Compare Predictions of fixed Re and mstar
%Data location from free case:
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1'}
mesh_plot={'Adapt_mode_Hmax10_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]};

Re=60;%20 21 31 35 40 45 60
mstar=20;

[Ustar_impedance,lambda_r_impedance,U_free,lambda_r_Free,U_freeFLUID,lambda_r_FreeFLUID]=SF_Impedance_Treatement(Re,mstar,folder_plot);

%Figures:
figure, hold on
plot(Ustar_impedance,lambda_r_impedance); %0.5 é por causa da differenca L e CL
plot(U_free,lambda_r_Free);
plot(U_freeFLUID,lambda_r_FreeFLUID);
plot(Ustar_impedance,Ustar_impedance*0,'k--');
title(['Comparision Impedance-Based-Method and Free Results for Re=' num2str(Re)]);
xlabel('U^*'); ylabel('\lambda_r');
legend('Impedance Based Predictions','Data from free Case');


%% OLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Re=40
% For m=100
Re=40;
mstar=1;

%Load Data Forced-Case:
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe_total=importFFdata(all_data_stored_file);
%Data Forced case:
Zr=2*real(dataRe_total.Lift);
Zi=2*imag(dataRe_total.Lift);
Omega_f=dataRe_total.OMEGAtab;

Ustar_impedance=sqrt((4*(pi^2))./(Omega_f.*(Omega_f+2./(pi.*mstar).*Zi)));
lambda_r_impedance=-Zr./(0.5*mstar*pi);

% Load data from Free case:
FreeCase=load([folder_plot{1} 'Re' num2str(Re) '/mstar' num2str(mstar) '/02modeSTRUCTURE_data.mat']);
lambda_r_Free=real(FreeCase.sigma_tab);
U_free=FreeCase.U_star;

%Figures:
figure, hold on
plot(Ustar_impedance,lambda_r_impedance); %0.5 é por causa da differenca L e CL
plot(U_free,lambda_r_Free);
plot([5 20],[0 0],'k--')
title('Comparision Impedance-Based-Method and Free Results for Re=40');
xlabel('U^*'); ylabel('\lambda_r or Z_r');
legend('Impedance Based Predictions','Data from free Case')





 %% Limit for mstar=50 as function of Re
 
 %estou a fazer isto...
 mstar=50;
 Re_tab=[20:2:45];
 curve_zone_unstable=[];
 
 
 for Re=Re_tab
     all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
     dataRe_total=importFFdata(all_data_stored_file);
     
     for i=1:((size(dataRe_total.OMEGA_tab,2))-1)
         if (real(dataRe_total.Lift(i))*real(dataRe_total.Lift(i+1))<0)
             curve_zone_unstable=[curve_zone_unstable [ Re dataRe_total.OMEGA_tab(i) imag(dataRe_total.Lift(i)) ]' ];
         end
     end
 end
 
 Re_f=curve_zone_unstable(1,:);
 Omega_f=dataRe_total.OMEGAtab(index);
 Zi_f=curve_zone_unstable(3,:);

 
 
 U_star_calculated=sqrt((4*pi^2)./(dataRe_total.OMEGA_tab.*(dataRe_total.OMEGA_tab+2./(pi.*mstar).*dataRe_total.Lift)));
 
 figure;
 scatter(curve_zone_unstable(1,:), U_star_calculated );