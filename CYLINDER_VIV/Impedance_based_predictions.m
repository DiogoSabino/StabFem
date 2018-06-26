%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Impedance method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Limit Rec1

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
DataFreeCase=load('./Impedance_Treatement/FreeCase.mat');
m_free=DataFreeCase.Point_TAB_Rec1(2,:);
U_free=DataFreeCase.Point_TAB_Rec1(3,:);

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


%% Courbe que David demande...
% For Re=21

Re=21;
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe_total=importFFdata(all_data_stored_file);

Zr=real(dataRe_total.Lift);
Zi=imag(dataRe_total.Lift);
Omega_f=dataRe_total.OMEGAtab;

mstar=100;
Ustar_impedance=sqrt((4*(pi^2))./(Omega_f.*(Omega_f+2./(pi.*mstar).*Zi)));

% Load data from Free case:
DataFreeCase=load('./Impedance_Treatement/FreeCaseRe21mstar100.mat');
lambda_r_Free=real(DataFreeCase.DataFree_curveUstar_Re21_m100.sigma_tab);
U_free=DataFreeCase.DataFree_curveUstar_Re21_m100.U_star;

figure, hold on
plot(Ustar_impedance,-Zr./(300*pi));
plot(U_free,lambda_r_Free);
plot([8.5 10.5],[0 0],'k--')
title('Comparision Impedance-Based-Method and Free Results for Re=21');
xlabel('U^*'); ylabel('\lambda_r or Z_r');
legend('Impedance Based Predictions','Data from free Case')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Re=40

Re=40;
all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe_total=importFFdata(all_data_stored_file);

Zr=real(dataRe_total.Lift);
Zi=imag(dataRe_total.Lift);
Omega_f=dataRe_total.OMEGAtab;

mstar=100;
Ustar_impedance=sqrt((4*(pi^2))./(Omega_f.*(Omega_f+2./(pi.*mstar).*Zi)));

% Load data from Free case:
DataFreeCase=load('./Impedance_Treatement/FreeCaseRe40mstar100.mat');
lambda_r_Free=real(DataFreeCase.DataFree_curveUstar_Re40_m100.sigma_tab);
U_free=DataFreeCase.DataFree_curveUstar_Re40_m100.U_star;

figure, hold on
plot(Ustar_impedance,-Zr./(300*pi));
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