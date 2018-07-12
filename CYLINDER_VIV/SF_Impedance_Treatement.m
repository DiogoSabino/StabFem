function [Ustar_impedance,lambda_r_impedance,U_free,lambda_r_Free,U_freeFLUID,lambda_r_FreeFLUID]=SF_Impedance_Treatement(Re,mstar,folder_plot,formulation)

global ffdataharmonicdir

%Load Data Forced-Case:
all_data_stored_file=[ffdataharmonicdir formulation 'Forced_Harmonic2D_Re' num2str(Re) 'TOTAL.ff2m'];
dataRe_total=importFFdata(all_data_stored_file);
diffZ=load([ffdataharmonicdir formulation 'Forced_Re' num2str(Re) '_diff_Lift_Coeff.mat']);
%Data Forced case:
Zr=2*real(dataRe_total.Lift);
Zi=2*imag(dataRe_total.Lift);
Omega_f=dataRe_total.OMEGAtab;
dZi=diffZ.dZi; dZr=diffZ.dZr;
dOMEGA=diffZ.dOMEGA;

%dZi=[dZi ;dZi(end)]; dZr=[dZr; dZr(end)]; %DAVID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dOMEGA=[dOMEGA(1) ; dOMEGA]; %DAVID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dZi=[dZi(1) ; dZi]; dZr=[dZr(1) ; dZr]; 
dZi=dZi./dOMEGA;
dZr=dZr./dOMEGA;

A=-Zr;
B=pi*mstar+dZr+Zi./Omega_f;
C=dZi;

Ustar_impedance=sqrt((4*(pi^2))./(Omega_f.*(Omega_f+2./(pi.*mstar).*Zi)))+imag(A./(B+C.*1i));%-A.*C./(B.^2+C.^2);
lambda_r_impedance=real(A./(B+C.*1i));%A.*B./(B.^2+C.^2);



%Load Data Free-Case:
FreeCase=load([folder_plot{1} 'Re' num2str(Re) '/mstar' num2str(mstar) '/02modeSTRUCTURE_data.mat']);
lambda_r_Free=real(FreeCase.sigma_tab);
U_free=FreeCase.U_star;



FLUID_path=[folder_plot{1} 'Re' num2str(Re) '/mstar' num2str(mstar) '/03modeFLUID_data.mat'];
if(exist(FLUID_path)~=2)
    lambda_r_FreeFLUID=[];
    U_freeFLUID=[];
else
    disp('EXISTS!!')
    FreeCase=load(FLUID_path);
    lambda_r_FreeFLUID=real(FreeCase.sigma_tab);
    U_freeFLUID=FreeCase.U_star;
end



end