%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Impedance method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ffdatadir ffdataharmonicdir verbosity
run('../SOURCES_MATLAB/SF_Start.m');
ffdatadir = './WORK/';
verbosity=100;
%% CHOOSE folder for saving data:
%General_data_dir='v1/General_data/';
General_data_dir='v1/General_data_refined/';
%General_data_dir='v1/Rec1/';
%General_data_dir='v1/Rec2/';
%General_data_dir='v1/Rec3/';
%General_data_dir='v1/Limit_St0/';
%General_data_dir='v1/Limit_St_Infty/';

domain_parameters=[-50 50 50];
%domain_parameters=[-100 100 100]; %for the case when St->0
domain_identity={[ num2str(domain_parameters(1)) '_' num2str(domain_parameters(2)) '_' num2str(domain_parameters(3)) '/']};

%mesh_identity={'Adapt_mode_Hmax2_InterError_0.02/'};
%mesh_identity={'Adapt_sensibility_Hmax2_InterError_0.02/'};
%mesh_identity={'Adapt_sensibility_Hmax1_InterError_0.02/'};
mesh_identity={'Adapt_S_Hmax1_InterError_0.02/'};

savedata_dir={[ General_data_dir domain_identity{1} mesh_identity{1}]};
% path of the saved data for the harmonic case cylinder (repeated in macros.edp)
ffdataharmonicdir={[ffdatadir 'DATA_Forced_Cylinder/' savedata_dir{1} ]};

% Create path if it does not exist
if(exist(ffdataharmonicdir{1})~=7&&exist(ffdataharmonicdir{1})~=5)
    mysystem(['mkdir -p ' ffdataharmonicdir{1}]); %I read in internet that the '-p'(stands for parent) not always work in every shell...
end %It's compulsory: Creation of the directory

formulation='R'; %Put the formulation used
%Normally, this is the general name of the file, so dont change it
filename=[formulation 'Forced_Harmonic2D_Re']; %name WITHOUT the Re

%% Order epsilon 0 Predictions at Re=19.95 
%Done

load([ffdataharmonicdir{1}  formulation 'DATA.mat'],'Re_c1_FINAL_REc1','Zimin_FINAL_REc1','Omegamin_FINAL_REc1');

%Variation of mstar:
mstar=[0:0.05:19.9 20:0.1:50 60:10:1000];
Ustar_impedance=sqrt((4*(pi^2))./(Omegamin_FINAL_REc1.*(Omegamin_FINAL_REc1+2./(pi.*mstar).*2*Zimin_FINAL_REc1)));
%Plot predictions:
figure
plot(mstar,Ustar_impedance);
figure
loglog(mstar,Ustar_impedance);
figure
semilogx(mstar,Ustar_impedance);
%%%%filename_latex=['./Latex_data/Free/Re20/bIMPEDANCEmesh50_50_50_Re19p95.txt'];%DONE in relatorio
%%%%for index=1:size(mstar,2)
%%%%    str_latex=['(' num2str(mstar(index)) ',' num2str(Ustar_impedance(index)) ')'];
%%%%    dlmwrite(filename_latex,str_latex,'delimiter', '','-append' )
%%%%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Pour enquanto uso o 19.95
%Data location from free case:
General_data_dir_folder='./Final_Results_v20/';
domain_plot={'-50_50_50/'};
mesh_plot={'Adapt_S_Hmax1_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]};

mstar_free=[0.1:0.1:1 2:1:10 20:10:50 100:100:500 1000];
TAB_FREE=[];
%%%%filename_latex=['./Latex_data/Free/Re20/bFREEmesh50_50_50_Re19p95.txt'];%DONE in relatorio

for mf=mstar_free
    FreeCase=load([folder_plot{1} 'Re' num2str(19.95) '/mstar' num2str(mf) '/02modeSTRUCTURE_data.mat']);
    [lambda_r_Free,index_max]=max(real(FreeCase.sigma_tab));
    U_max=FreeCase.U_star(index_max);
    TAB_FREE=[TAB_FREE, [mf  U_max]'];
    %%%%str_latex=['(' num2str(mf) ',' num2str(U_max) ')'];
    %%%%dlmwrite(filename_latex,str_latex,'delimiter', '','-append' )
end


%Figure
figure;
plot(mstar,Ustar_impedance); hold on;
scatter(TAB_FREE(1,:),TAB_FREE(2,:),'r');

title('Predictions of Impedance Based Method');
xlabel('m^*'); ylabel('U^*');
legend('Impedance Based Predictions','Data from free Case')

%Figure with x-axis in log scale
figure;
semilogx(mstar,Ustar_impedance); hold on;
scatter(TAB_FREE(1,:),TAB_FREE(2,:),'r');

title('Predictions of Impedance Based Method (in log scale)');
xlabel('m^*'); ylabel('U^*');
legend('Impedance Based Predictions','Data from free Case')


%% Order 2 predictions

%Data location from free case:
General_data_dir_folder='./Final_Results_v20/';    %General_data_dir; % e.g.: './FOLDER_TOTO/'
domain_plot={'-50_50_50/'}; %domain_identity;        %e.g.:{'totodir1'}
mesh_plot={'Adapt_D_Hmax10_InterError_0.02/'};
folder_plot={[General_data_dir_folder  domain_plot{1} mesh_plot{1} ]};

Re=40;%20 21 31 35 40 45 60
mstar=100;

[Ustar_impedance,lambda_r_impedance,U_free,lambda_r_Free,U_freeFLUID,lambda_r_FreeFLUID]=SF_Impedance_Treatement(Re,mstar,folder_plot,filename);

%Figures:
figure, hold on
plot(Ustar_impedance,lambda_r_impedance);
plot(U_free,lambda_r_Free);
plot(U_freeFLUID,lambda_r_FreeFLUID);
plot(Ustar_impedance,Ustar_impedance*0,'k--');
title(['Comparision Impedance-Based-Method and Free Results for Re=' num2str(Re)]);
xlabel('U^*'); ylabel('\lambda_r');
legend('Impedance Based Predictions','Data from free Case');




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