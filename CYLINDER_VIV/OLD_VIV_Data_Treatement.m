%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is part of StabFem Project.
%   (set of programs in Freefem to perform global stability calculations  
%   in hydrodynamic stability)
%
%	File: Main_VIV_Data_Treatement.m
%   Contributours: David Fabre, Diogo Sabino ...
%   Last Modification: Diogo Sabino, 16 May 2018
%   
%	This script treats data already generated of flow aroung a rigid circular 
%   cilinder mounted in a spring-damped system to simulate VIV.
%
%   To give a good use of this script, run each section at a time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CE FICHIER CEST LE BORDEL, IL FAUT QUE JE M'EN OCUPE

%% Forced case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%à faire proprement

%% Free Movement 1 transverse DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOMAINE VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RE60 mstar20
clear all
close all

%import data from NAVROSE and plot it
Navrose_Re60_mstar20 = importdata('./Navrose_Data/RE60_M20_real.csv');

figure; hold on
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,2),'--*');
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,3),'--*','HandleVisibility','off');

%Plot our data
m50x50x50_FLUID=load('./Final_results_vm50x50x50/Re60/mstar20/03modeFLUI_spectrum.mat','sigma_tab','U_star');
%m50x50x50_STRU=load('./Final_results_vm50x50x50/Re60/mstar20/03modeSTRU_spectrum.mat','sigma_tab','U_star','HandleVisibility','off');
plot(m50x50x50_FLUID.U_star,real(m50x50x50_FLUID.sigma_tab),'-*r');
%plot(m50x50x50_STRU.U_star,real(m50x50x50_STRU.sigma_tab),'-*r');
m40x80x40_FLUID=load('./Final_results_vm40x80x40/Re60/mstar20/03modeFLUI_spectrum.mat','sigma_tab','U_star');
%m40x80x40_STRU=load('./Final_results_vm40x80x40/Re60/mstar20/03modeSTRU_spectrum.mat','sigma_tab','U_star','HandleVisibility','off');
plot(m40x80x40_FLUID.U_star,real(m40x80x40_FLUID.sigma_tab),'-*b');
%plot(m40x80x40_STRU.U_star,real(m40x80x40_STRU.sigma_tab),'-*b');
m80x160x80_FLUID=load('./Final_results_vm80x160x80/Re60/mstar20/03modeFLUI_spectrum.mat','sigma_tab','U_star');
%m80x160x80_STRU=load('./Final_results_vm80x160x80/Re60/mstar20/03modeSTRU_spectrum.mat','sigma_tab','U_star','HandleVisibility','off');
plot(m80x160x80_FLUID.U_star,real(m80x160x80_FLUID.sigma_tab),'-*m');
%plot(m80x160x80_STRU.U_star,real(m80x160x80_STRU.sigma_tab),'-*m');
%m80x160x80_split_FLUID=load('./Final_results_vm80x160x80_split/Re60/mstar20/03modeFLUI_spectrum.mat','sigma_tab','U_star');
%m80x160x80_split_STRU=load('./Final_results_vm80x160x80_split/Re60/mstar20/03modeSTRU_spectrum.mat','sigma_tab','U_star','HandleVisibility','off');
%plot(m80x160x80_split_FLUID.U_star,real(m80x160x80_split_FLUID.sigma_tab),'-*y');
%plot(m80x160x80_split_STRU.U_star,real(m80x160x80_split_STRU.sigma_tab),'-*y');

title('Amplification Rate: Domaine variation for Re=60, mstar=20, Hmax=10, InterpError=0.02')
legend('Navrose Data D=[-50:50]x[0:50]','D=[-50:50]x[0:50]','D=[-40:80]x[0:40]','D=[-80:160]x[0:80]')
%saveas(gcf,'./General_images/domain_variationRE60M20.png');
%End RE60M20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RE40
%clear all
close all

%import data from NAVROSE and plot it
Navrose_Re40_mstar10 = importdata('./Navrose_Data/RE40_M10_real.csv');

figure; hold on
plot(Navrose_Re40_mstar10.data(:,1),Navrose_Re40_mstar10.data(:,2),'--*');
%plot(Navrose_Re40_mstar10.data(:,1),Navrose_Re40_mstar10.data(:,3),'--*');

%Plot our data
m50x50x50_FLUID=load('./Final_results_vm50x50x50/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
plot(m50x50x50_FLUID.U_star,real(m50x50x50_FLUID.sigma_tab),'-*r');
m40x80x40_FLUID=load('./Final_results_vm40x80x40/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
plot(m40x80x40_FLUID.U_star,real(m40x80x40_FLUID.sigma_tab),'-*b');
m80x160x80_FLUID=load('./Final_results_vm80x160x80/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
plot(m80x160x80_FLUID.U_star,real(m80x160x80_FLUID.sigma_tab),'-*m');
m40x80x40_mesh_different_FLUID=load('./Final_results_vm40x80x40_mesh_different/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
plot(m40x80x40_mesh_different_FLUID.U_star,real(m40x80x40_mesh_different_FLUID.sigma_tab),'-*y');
m40x80x40_mesh_different_2_FLUID=load('./Final_results_vm40x80x40_mesh_different_2/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
plot(m40x80x40_mesh_different_2_FLUID.U_star,real(m40x80x40_mesh_different_2_FLUID.sigma_tab),'-*c');


title('Amplification Rate: Domaine variation for Re=40, mstar=10, Hmax=10, InterpError=0.02')
legend('Navrose Data D=[-50:50]x[0:50]','D=[-50:50]x[0:50]','D=[-40:80]x[0:40]','D=[-80:160]x[0:80]','D=[-40:80]x[0:40], finner mesh','D=[-40:80]x[0:40], w/out mode adaptation')
%saveas(gcf,'./General_images/domain_variationRE40M10.png');
%End RE40


%% Comparasion with NAVROSE
clear all
close all
%RE60 m20
open('./Final_results_v13/Re60/mstar20/04mode_both.fig');

Navrose_Re60_mstar20 = importdata('./Navrose_Data/RE60_M20_real.csv');
subplot(2,1,1); hold on;
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,2),'--*');
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,3),'--*');
legend('StrucMODE','FluiMODE','FluiMODE Navrose','StrucMODE Navrose')


Navrose_Re60_mstar20 = importdata('./Navrose_Data/RE60_M20_imag.csv');
subplot(2,1,2); hold on;
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,2),'--*');
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,3),'--*');
legend('StrucMODE','FluiMODE','FluiMODE Navrose','StrucMODE Navrose')

%saveas(gcf,'./General_images/RE60_mstar20_comparasion.png');
%RE60 m5
open('./Final_results_v13/Re60/mstar5/04mode_both.fig');

Navrose_Re60_mstar5 = importdata('./Navrose_Data/RE60_M5_real.csv');
subplot(2,1,1); hold on;
plot(Navrose_Re60_mstar5.data(:,1),Navrose_Re60_mstar5.data(:,2),'--*');
plot(Navrose_Re60_mstar5.data(:,1),Navrose_Re60_mstar5.data(:,3),'--*');
legend('FEM2MODE','FEM1MODE','FEM1MODE Navrose','FEM2MODE Navrose')

Navrose_Re60_mstar5 = importdata('./Navrose_Data/RE60_M5_imag.csv');
subplot(2,1,2); hold on;
plot(Navrose_Re60_mstar5.data(:,1),Navrose_Re60_mstar5.data(:,2),'--*');
plot(Navrose_Re60_mstar5.data(:,1),Navrose_Re60_mstar5.data(:,3),'--*');
legend('FEM2MODE','FEM1MODE','FEM1MODE Navrose','FEM2MODE Navrose')

%saveas(gcf,'./General_images/RE60_mstar5_comparasion.png');

clear all
close all

%% Analysis for different masses
%RE40
%open('./Final_results_v13/Re40/mstar10/06_Navrose_p578.fig');

mstar100_FLUID=load('./Final_results_v13/Re40/mstar100/03modeFLUI_spectrum.mat','sigma_tab','U_star');
mstar20_FLUID=load('./Final_results_v13/Re40/mstar20/03modeFLUI_spectrum.mat','sigma_tab','U_star');
mstar15_FLUID=load('./Final_results_v13/Re40/mstar15/03modeFLUI_spectrum.mat','sigma_tab','U_star');
mstar10_FLUID=load('./Final_results_v13/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');

mstar100_STRU=load('./Final_results_v13/Re40/mstar100/02modeSTRU_spectrum.mat','sigma_tab','U_star');
mstar20_STRU=load('./Final_results_v13/Re40/mstar20/02modeSTRU_spectrum.mat','sigma_tab','U_star');
mstar15_STRU=load('./Final_results_v13/Re40/mstar15/02modeSTRU_spectrum.mat','sigma_tab','U_star');
mstar10_STRU=load('./Final_results_v13/Re40/mstar10/02modeSTRU_spectrum.mat','sigma_tab','U_star');

Navrose_Re40_mstar10 = importdata('./Navrose_Data/RE40_M10_real.csv');



figure; hold on
%subplot(2,1,1); hold on
plot(mstar100_FLUID.U_star,real(mstar100_FLUID.sigma_tab),'-*r');
plot(mstar20_FLUID.U_star,real(mstar20_FLUID.sigma_tab),'-*b');
plot(mstar15_FLUID.U_star,real(mstar15_FLUID.sigma_tab),'-*m');
plot(mstar10_FLUID.U_star,real(mstar10_FLUID.sigma_tab),'-*g');
plot(mstar100_STRU.U_star,real(mstar100_STRU.sigma_tab),'-or','HandleVisibility','off','MarkerSize',2);
plot(mstar20_STRU.U_star,real(mstar20_STRU.sigma_tab),'-ob','HandleVisibility','off','MarkerSize',2);
plot(mstar15_STRU.U_star,real(mstar15_STRU.sigma_tab),'-om','HandleVisibility','off','MarkerSize',2);
plot(mstar10_STRU.U_star,real(mstar10_STRU.sigma_tab),'-og','HandleVisibility','off','MarkerSize',2);
title('Amplification Rate: Mass variation for Re=40 (*: Mode FLUID ; °: Mode STRUCTURE)')
legend('mstar=100','mstar=20','mstar=15','mstar=10')

plot(Navrose_Re40_mstar10.data(:,1),Navrose_Re40_mstar10.data(:,2),'--*');
plot(Navrose_Re40_mstar10.data(:,1),Navrose_Re40_mstar10.data(:,3),'--*');
% subplot(2,1,2); hold on
% plot(mstar100_FLUID.U_star,imag(mstar100_FLUID.sigma_tab).*mstar100_FLUID.U_star/2/pi,'-r');
% plot(mstar20_FLUID.U_star,imag(mstar20_FLUID.sigma_tab).*mstar20_FLUID.U_star/2/pi,'-b');
% plot(mstar15_FLUID.U_star,imag(mstar15_FLUID.sigma_tab).*mstar15_FLUID.U_star/2/pi,'-m');
% plot(mstar10_FLUID.U_star,imag(mstar10_FLUID.sigma_tab).*mstar10_FLUID.U_star/2/pi,'-g');
% plot(mstar100_STRU.U_star,imag(mstar100_STRU.sigma_tab).*mstar100_FLUID.U_star/2/pi,'-r','HandleVisibility','off');
% plot(mstar20_STRU.U_star,imag(mstar20_STRU.sigma_tab).*mstar20_FLUID.U_star/2/pi,'-b','HandleVisibility','off');
% plot(mstar15_STRU.U_star,imag(mstar15_STRU.sigma_tab).*mstar15_FLUID.U_star/2/pi,'-m','HandleVisibility','off');
% plot(mstar10_STRU.U_star,imag(mstar10_STRU.sigma_tab).*mstar10_FLUID.U_star/2/pi,'-g','HandleVisibility','off');

%% à la main Comparasion with NAVROSE D33
%clear all
%close all
%RE60 m20
open('./Final_results_tests_m_different/Re60/mstar20/03modeFLUI.fig');
%é mais metodico de fazer load do ficheiro mat, do que abrir a figura guardada

mstar20_FLUID=load('./Final_results_v13/Re60/mstar20/03modeFLUI_spectrum.mat','sigma_tab','U_star');
subplot(2,1,1); hold on;
plot(mstar20_FLUID.U_star,real(mstar20_FLUID.sigma_tab));


Navrose_Re60_mstar20 = importdata('./Navrose_Data/RE60_M20_real.csv');
subplot(2,1,1); hold on;
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,2),'--*');
%legend('FluiMODE with adapt mesh each 10iterations','FluiMODE just 1 initial adapt mesh','FluiMODE Navrose')

%% à la main Comparasion with NAVROSE  D33


%Load data to compare
vm40x80x40_Re40_mstar10_FLUID=load('./Final_results_vm40x80x40/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
vm50x50x50_Re40_mstar10_FLUID=load('./Final_results_vm50x50x50/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
vm40x80x40ADAPT_10_Re40_mstar10_FLUID=load('./Final_results_vm40x80x40ADAPT_10/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
vm40x80x40ADAPT_4_Re40_mstar10_FLUID=load('./Final_results_vm40x80x40ADAPT_4/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');
vm50x50x50_ADAPT_4_Re40_mstar10_FLUID=load('./Final_results_vm50x50x50_ADAPT_4/Re40/mstar10/03modeFLUI_spectrum.mat','sigma_tab','U_star');

Navrose_Re40_mstar10 = importdata('./Navrose_Data/RE40_M10_real.csv'); %meter os bons nomes

%Plot the desired data
figure; hold on

plot(vm40x80x40_Re40_mstar10_FLUID.U_star,real(vm40x80x40_Re40_mstar10_FLUID.sigma_tab));
plot(vm50x50x50_Re40_mstar10_FLUID.U_star,real(vm50x50x50_Re40_mstar10_FLUID.sigma_tab));
plot(vm40x80x40ADAPT_10_Re40_mstar10_FLUID.U_star,real(vm40x80x40ADAPT_10_Re40_mstar10_FLUID.sigma_tab));
plot(vm40x80x40ADAPT_4_Re40_mstar10_FLUID.U_star,real(vm40x80x40ADAPT_4_Re40_mstar10_FLUID.sigma_tab));
plot(vm50x50x50_ADAPT_4_Re40_mstar10_FLUID.U_star,real(vm50x50x50_ADAPT_4_Re40_mstar10_FLUID.sigma_tab));

%Navrose data
plot(Navrose_Re40_mstar10.data(:,1),Navrose_Re40_mstar10.data(:,2),'--*'); %ver qual é que é o certo

%title('Amplification Rate: Mass variation for Re=40 (*: Mode FLUID ; °: Mode STRUCTURE)')
legend('without ADAPTMESH40x80x40','without ADAPTMESH50x50x50','with ADAPT40x80x40 at each 10it','with ADAPT40x80x40 at each 4it','with ADAPT50x50x50 at each 4it','Navrose')

%% Day when I understand the terms associated to C_F

v20_Re60_mstar20_FLUID=load('./Final_results_v13/Re60/mstar20/03modeFLUI_spectrum.mat','sigma_tab','U_star');
v21_Re60_mstar20_FLUID=load('./???????/Re60/mstar20/03modeFLUID_data.mat','sigma_tab','U_star');
Navrose_Re60_mstar20 = importdata('./Navrose_Data/RE60_M20_real.csv'); %meter os bons nomes

figure; hold on

plot(v20_Re60_mstar20_FLUID.U_star,real(v20_Re60_mstar20_FLUID.sigma_tab));
plot(v21_Re60_mstar20_FLUID.U_star,real(v21_Re60_mstar20_FLUID.sigma_tab));
%Navrose data
plot(Navrose_Re60_mstar20.data(:,1),Navrose_Re60_mstar20.data(:,2),'--*'); %ver qual é que é o certo

legend('Dabre et Diogo','essay','Navrose');


