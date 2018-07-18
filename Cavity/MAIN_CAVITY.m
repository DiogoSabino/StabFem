clear all;
close all;
run('../SOURCES_MATLAB/SF_Start.m');
figureformat='png'; AspectRatio = 0.56; % for figures
   
disp(' STARTING ADAPTMESH PROCEDURE : ');    
disp(' ');
disp(' ');    
verbosity=100;
bf=SF_Init('Mesh_Cavity.edp'); %Do it with input parameters



bf=SF_BaseFlow(bf,'Re',500);
bf=SF_BaseFlow(bf,'Re',1000);
bf=SF_BaseFlow(bf,'Re',1500);
bf=SF_BaseFlow(bf,'Re',3000);
bf=SF_BaseFlow(bf,'Re',4140);

bf=SF_Adapt(bf,'Hmax',0.25);

bf=SF_BaseFlow(bf,'Re',4140); % not needed, I think
disp(' ');
disp('mesh adaptation to SENSIBILITY : ')

[ev,em] = SF_Stability(bf,'shift',0.0+7.5i,'nev',1,'type','S');

[bf,em]=SF_Adapt(bf,em,'Hmax',0.25);
%% HB_O1 for Mode wiht Omega=10.5 at Re=4570 (2nd branch)

bf=SF_BaseFlow(bf,'Re',4570);
[ev,em] = SF_Stability(bf,'shift',0.0+10.5i,'nev',1,'type','D');
[meanflow,mode] = SF_SelfConsistentDirect(bf,em,'sigma',0.,'Re',4570); 
omega_HB1 = [imag(mode.lambda)]; Aenergy_HB1  = [mode.AEnergy];
system(['cp ' ffdatadir  'SelfConsistentMode.txt ' './Results_v1/Guess_for_2nd_branchHB1/Re4570SelfConsistentMode.txt' ]);
system(['cp ' ffdatadir  'MeanFlow.txt ' './Results_v1/Guess_for_2nd_branchHB1/Re4570MeanFlow.txt' ]);

Re_tab=[4600:100:7000];
for Re=Re_tab
[meanflow,mode] = SF_SelfConsistentDirect(meanflow,mode,'Re',Re);
omega_HB1 = [omega_HB1 imag(mode.lambda)];
Aenergy_HB1  = [Aenergy_HB1 mode.AEnergy];
system(['cp ' ffdatadir  'MeanFlow.txt ' './Results_v1/Guess_for_2nd_branchHB1/Re' num2str(Re) 'MeanFlow.txt' ]);
system(['cp ' ffdatadir  'SelfConsistentMode.txt ' './Results_v1/Guess_for_2nd_branchHB1/Re' num2str(Re) 'SelfConsistentMode.txt' ]);
end

if(1==0)
    %Re_tab=[4570 Re_tab];
    %save('./Results_v1/Guess_for_2nd_branchHB1/data.mat','Aenergy_HB1','omega_HB1','Re_tab')
    figure; hold on;
    subplot(1,2,1)
    plot(Re_tab,Aenergy_HB1)
        subplot(1,2,2)
    plot(Re_tab,omega_HB1)

end

%% HB_O2 for Mode wiht Omega=7.5 at Re=4160 (1st branch)
bf=SF_BaseFlow(bf,'Re',4160);
[ev,em] = SF_Stability(bf,'shift',0.0+7.5i,'nev',1,'type','D');
[meanflow,mode] = SF_SelfConsistentDirect(bf,em,'sigma',0.,'Re',4160);
system(['cp ' ffdatadir 'SelfConsistentMode_guess.txt '  ffdatadir 'SecondHarmonicMode_guess.txt ']);

[meanflow,mode,mode2] = SF_HarmonicBalance_Ordre2(meanflow,mode,mode);
system(['cp ' ffdatadir  'MeanFlow.txt ' './Results_v1/Guess_for_1nd_branchHB2/Re4570MeanFlow.txt' ]);
system(['cp ' ffdatadir  'SelfConsistentMode.txt ' './Results_v1/Guess_for_1nd_branchHB2/Re4570SelfConsistentMode.txt' ]);
system(['cp ' ffdatadir  'SecondHarmonicMode.txt ' './Results_v1/Guess_for_1nd_branchHB2/Re4570SecondHarmonicMode.txt' ]);

omega_HB2 = [imag(mode.lambda)]; Aenergy_HB2  = [sqrt(mode.AEnergy^2+ mode2.AEnergy^2)];

Re_tab=[4200:100:7000];
for Re=Re_tab
    [meanflow,mode,mode2] = SF_HarmonicBalance_Ordre2(meanflow,mode,mode2,'Re',Re);
    omega_HB2 = [omega_HB2 imag(mode.lambda)];
    Aenergy_HB2 = [Aenergy_HB2 sqrt(mode.AEnergy^2+ mode2.AEnergy^2)];
    system(['cp ' ffdatadir  'MeanFlow.txt ' './Results_v1/Guess_for_1nd_branchHB2/Re' num2str(Re) 'MeanFlow.txt' ]);
    system(['cp ' ffdatadir  'SelfConsistentMode.txt ' './Results_v1/Guess_for_1nd_branchHB2/Re' num2str(Re) 'SelfConsistentMode.txt' ]);
    system(['cp ' ffdatadir  'SecondHarmonicMode.txt ' './Results_v1/Guess_for_1nd_branchHB2/Re' num2str(Re) 'SecondHarmonicMode.txt' ]);
end
%Re_tab=[4160 Re_tab];
%save('./Results_v1/Guess_for_1nd_branchHB2/data.mat','omega_HB2','Aenergy_HB2','Re_tab' )

%% HB_O2 for Mode wiht Omega=10.5 at Re=4570 (2nd branch)

bf=SF_BaseFlow(bf,'Re',4570);
[ev,em] = SF_Stability(bf,'shift',0.0+10.5i,'nev',1,'type','D');
[meanflow,mode] = SF_SelfConsistentDirect(bf,em,'sigma',0.,'Re',4570);

system(['cp ' ffdatadir 'SelfConsistentMode_guess.txt '  ffdatadir 'SecondHarmonicMode_guess.txt ']);

%[meanflow,mode,mode2] = SF_HarmonicBalance_Ordre2(meanflow1,mode1,mode1);
[meanflow,mode,mode2] = SF_HarmonicBalance_Ordre2(meanflow,mode,mode);
system(['cp ' ffdatadir  'MeanFlow.txt ' './Results_v1/Guess_for_2nd_branchHB2/Re4570MeanFlow.txt' ]);
system(['cp ' ffdatadir  'SelfConsistentMode.txt ' './Results_v1/Guess_for_2nd_branchHB2/Re4570SelfConsistentMode.txt' ]);
system(['cp ' ffdatadir  'SecondHarmonicMode.txt ' './Results_v1/Guess_for_2nd_branchHB2/Re4750SecondHarmonicMode.txt' ]); %%%ERROOO

omega_HB2 = [imag(mode.lambda)]; Aenergy_HB2  = [sqrt(mode.AEnergy^2+ mode2.AEnergy^2)];


%Re_tab=[4600:50:7000];
%Re_tab=[4700 4900 5100 5200 5400 5600 5700 5800:100:7000];
%Re_tab=[6200:100:7000];
Re_tab=[4500:-100:4200 4160];
for Re=Re_tab
    [meanflow,mode,mode2] = SF_HarmonicBalance_Ordre2(meanflow,mode,mode2,'Re',Re);
    omega_HB2 = [omega_HB2 imag(mode.lambda)];
    Aenergy_HB2 = [Aenergy_HB2 sqrt(mode.AEnergy^2+ mode2.AEnergy^2)];
    system(['cp ' ffdatadir  'MeanFlow.txt ' './Results_v1/Guess_for_2nd_branchHB2/Re' num2str(Re) 'MeanFlow.txt' ]);
    system(['cp ' ffdatadir  'SelfConsistentMode.txt ' './Results_v1/Guess_for_2nd_branchHB2/Re' num2str(Re) 'SelfConsistentMode.txt' ]);
    system(['cp ' ffdatadir  'SecondHarmonicMode.txt ' './Results_v1/Guess_for_2nd_branchHB2/Re' num2str(Re) 'SecondHarmonicMode.txt' ]);
end

%% When something goes wrong: to retreive the files...
if(1==0)
    Re=4570;
    system(['cp ./Results_v1/Guess_for_2nd_branchHB2/Re' num2str(Re) 'MeanFlow.txt ' ffdatadir  'MeanFlow.txt ' ]);
    system(['cp ./Results_v1/Guess_for_2nd_branchHB2/Re' num2str(Re) 'SelfConsistentMode.txt ' ffdatadir  'SelfConsistentMode.txt ' ]);
    system(['cp ./Results_v1/Guess_for_2nd_branchHB2/Re' num2str(Re) 'SecondHarmonicMode.txt ' ffdatadir  'SecondHarmonicMode.txt ' ]);
end
%save('./Results_v1/Guess_for_2nd_branch/secondpart_from6200.mat','omega_HB2','Aenergy_HB2','Re_tab' )

%% Data post-traetement

%image Amplitude 2nd Branch
%image omega   2nd Branch

Meliga_E_SC = importdata('./Litterature_data/DATA/MeligaRe_Energy_SC.csv');
Meliga_E_DNS = importdata('./Litterature_data/DATA/MeligaRe_Energy_DNS.csv');
Meliga_Omega_SC = importdata('./Litterature_data/DATA/MeligaRe_Omega_SC.csv');
Meliga_Omega_DNS = importdata('./Litterature_data/DATA/MeligaRe_Omega_DNS.csv');
Diogo_HB1_Br2=load('./Results_v1/Guess_for_2nd_branchHB1/data.mat');
Diogo_HB2_Br2=load('./Results_v1/Guess_for_2nd_branchHB2/data.mat');

figure; hold on;
subplot(1,2,1);hold on;
plot(Meliga_E_SC.data(:,1),Meliga_E_SC.data(:,2),'ko','LineWidth',2) %SC MELIGA
Meliga_DNS_smooth_x = Meliga_E_DNS.data(1,1):0.25:Meliga_E_DNS.data(end,1);
Meliga_DNS_smooth_y = interp1(Meliga_E_DNS.data(:,1),Meliga_E_DNS.data(:,2),Meliga_DNS_smooth_x);
plot(Meliga_DNS_smooth_x,Meliga_DNS_smooth_y,'k-.','LineWidth',2) %DNS MELIGA

plot(Diogo_HB1_Br2.Re_tab,Diogo_HB1_Br2.Aenergy_HB1./sqrt(2),'r--','LineWidth',2) %DIOGO
plot(Diogo_HB2_Br2.Re_tab(1,4:end),Diogo_HB2_Br2.Aenergy_HB2(1,4:end)./sqrt(2),'b-','LineWidth',2) %DIOGO
%Details:
text(4000,0.007,'M17: SC','FontSize', 15); %DIOGO ICI
text(6400,0.039,'M17: DNS','FontSize', 15); %DIOGO ICI
text(6450,0.055,'HB1','FontSize', 15); %DIOGO ICI
text(6550,0.0455,'HB2','FontSize', 15); %DIOGO ICI
text(4050,0.056,'(a)','FontSize', 20); %DIOGO ICI
xlim([4000 7050]); %ylim([-0.4 0.8])
set(gca,'fontsize',15); box on
%legend('P. Meliga 2017: SC','P. Meliga 2017: DNS','Present HB1','Present HB2','Location','southeast')
xlabel('Re','FontSize', 20);ylabel('A','FontSize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2);hold on;
plot(Meliga_Omega_SC.data(:,1),Meliga_Omega_SC.data(:,2),'ko','LineWidth',2)
Meliga_DNS_smooth_x = Meliga_Omega_DNS.data(1,1):0.25:Meliga_Omega_DNS.data(end,1);
Meliga_DNS_smooth_y = interp1(Meliga_Omega_DNS.data(:,1),Meliga_Omega_DNS.data(:,2),Meliga_DNS_smooth_x);
plot(Meliga_DNS_smooth_x,Meliga_DNS_smooth_y,'k-.','LineWidth',2)

plot(Diogo_HB1_Br2.Re_tab,Diogo_HB1_Br2.omega_HB1,'r--','LineWidth',2) %DIOGO hb1
plot(Diogo_HB2_Br2.Re_tab(1,4:end),Diogo_HB2_Br2.omega_HB2(1,4:end),'b-','LineWidth',2) %DIOGO 

%Details:
%legend('P. Meliga 2017: SC','P. Meliga 2017: DNS','Present HB1','Present HB2','Location','southeast')
xlabel('Re','FontSize', 25);ylabel('\fontsize{25} \omega','FontSize',25);
text(4050,10.43,'M17: SC','FontSize', 15); %DIOGO ICI
text(5770,11.30,'M17: DNS','FontSize', 13); %DIOGO ICI
text(6550,11.69,'HB1','FontSize', 15); %DIOGO ICI
text(6450,11.35,'HB2','FontSize', 15); %DIOGO ICI
text(4050,11.7,'(b)','FontSize', 20); %DIOGO ICI
xlim([4000 7050]); box on; x0=10; y0=10; width=1800; height=400;
set(gca,'fontsize',15); set(gcf,'units','points','position',[x0,y0,width,height])
saveas(gcf,'Branch2_Energy_Omega.eps','epsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image omega: 1b of Meliga with the 2 branches

figure; hold on;
Diogo_HB2_Br1=load('./Results_v1/Guess_for_1nd_branchHB2/data.mat');
Meliga_Om_DNS_2Br = importdata('./Litterature_data/DATA/MeligaRe_Omega_DNS_2Br.csv');

plot(Diogo_HB2_Br2.Re_tab(1,1:4),Diogo_HB2_Br2.omega_HB2(1,1:4),'r--','HandleVisibility','off') %DIOGO BRANCH 2
plot(Diogo_HB2_Br2.Re_tab(1,4:end),Diogo_HB2_Br2.omega_HB2(1,4:end),'r-','LineWidth',2) %DIOGO BRANCH 2
plot(Diogo_HB2_Br1.Re_tab(1,1:4),Diogo_HB2_Br1.omega_HB2(1,1:4),'b-','LineWidth',2) %DIOGO BRANCH 1
plot(Diogo_HB2_Br1.Re_tab(1,4:end),Diogo_HB2_Br1.omega_HB2(1,4:end),'b--','HandleVisibility','off') %DIOGO BRANCH 1
plot(Meliga_Om_DNS_2Br.data(:,1),Meliga_Om_DNS_2Br.data(:,2),'ko','LineWidth',2) %2Br MELIGA


%legend('Present Br2','Present Br1','P. Meliga 2017: DNS','Location','southeast')
plot([4400 4400],[7 12],'k--','HandleVisibility','off')
%Details:
xlabel('Re','FontSize', 20); ylabel('\fontsize{25} \omega','FontSize', 25);
box on
text(4810,10.85,'M17: DNS','FontSize', 15); %DIOGO ICI
text(4100,10.5,'Branch 2','FontSize', 15); %DIOGO ICI
text(6000,8.2,'Branch 1','FontSize', 15); %DIOGO ICI
x0=10;
y0=410;
width=950;
height=400;
set(gca,'fontsize',15) 
set(gcf,'units','points','position',[x0,y0,width,height])
saveas(gcf,'Omega_Branches_1_2.eps','epsc')
%% Trash diogo

Meliga_E_SC = importdata('./Litterature_data/DATA/MeligaRe_Energy_SC.csv');
figure; hold on
plot(Meliga_E_SC.data(:,1),Meliga_E_SC.data(:,2),'o')

%FIRSTPART=load('./Results_v1/Resultsforthe_second_branch_firstpart.mat');
%plot([4570 FIRSTPART.Re_tab(1,1:25)],FIRSTPART.Aenergy_HB2(1,1:26)./sqrt(2))

%Re_tab=[4570 FIRSTPART.Re_tab(1,1:25)];
%Aenergy_HB2=FIRSTPART.Aenergy_HB2(1,1:26);
%omega_HB2=FIRSTPART.omega_HB2(1,1:26);

%save('./Results_v1/Guess_for_2nd_branchHB2/data_fromRe4570_toRe5800.mat','Re_tab','Aenergy_HB2','omega_HB2')

%plot(SECONDPART.Re_tab(1,:),SECONDPART.Aenergy_HB2(1,end-8:end)./sqrt(2))
%Re_tab=SECONDPART.Re_tab(1,:);
%Aenergy_HB2=SECONDPART.Aenergy_HB2(1,end-8:end);
%omega_HB2=SECONDPART.Aenergy_HB2(1,end-8:end);
 
%save('./Results_v1/Guess_for_2nd_branchHB2/data_fromRe6200_toRe7000.mat','Re_tab','Aenergy_HB2','omega_HB2')
 
%save('./Results_v1/Guess_for_2nd_branchHB2/data_fromRe5900_toRe6100.mat','pequeroRe','peauenoAenergia','peauenoomega')
 
%save('./Results_v1/Guess_for_2nd_branchHB2/data.mat','Re_tab','Aenergy_HB2','omega_HB2')
 





