%  Instability of the wake of a cylinder with STABFEM  
%
%  this script demonstrates the main functionalities of StabFem 
%  1/ Generation of an adapted mesh and base flow
%  2/ Eigenvalue computation
%  3/ New mesh adaptation to eigenmode

clear all; close all;
%set the global variables needed by the drivers
run('../SOURCES_MATLAB/SF_Start.m');


disp(' STARTING ADAPTMESH PROCEDURE : ');    
disp(' ');
disp(' LARGE MESH : [-40:80]x[0:40] ');
disp(' ');    
bf=SF_Init('Mesh_Cylinder.edp',[-40 80 40]);
plotFF(bf,'mesh');
pause(0.01);

bf=SF_BaseFlow(bf,'Re',1);
bf=SF_BaseFlow(bf,'Re',10);
bf=SF_BaseFlow(bf,'Re',60);
bf=SF_Adapt(bf);

bf
plotFF(bf,'ux');
pause(0.01);

disp(' ');

disp('Eigenvalue computation : compute 10 modes and plot the spectrum')
[ev,em] = SF_Stability(bf,'shift',0.04+0.76i,'nev',10,'type','D');
ev

disp('Eigenvalue computation : compute 1 mode, including structural sensitivity')
[ev,em] = SF_Stability(bf,'shift',0.04+0.76i,'nev',1,'type','S');


disp('mesh adaptation to SENSITIVITY and new computation of eigenvalues on improved mesh : ')
bf = SF_Adapt(bf,em);
[ev,em] = SF_Stability(bf,'shift',0.04+0.76i,'nev',1,'type','S');
ev
plotFF(em,'ux1');
plotFF(em,'sensitivity');
pause(0.01);

[ev,em] = SF_Stability(bf,'shift',0.04+0.76i,'nev',1,'type','D','plotspectrum','yes');



%%% DETERMINATION OF THE INSTABILITY THRESHOLD

bf=SF_BaseFlow(bf,'Re',50);
[ev,em] = SF_Stability(bf,'shift',+.75i,'nev',1,'type','D');
[bf,em]=SF_FindThreshold(bf,em);
Rec = bf.Re  
Cxc = bf.Cx; 
Lxc=bf.Lx;    Omegac=imag(em.lambda);

% Rec should be 47.6

%%% SelfConsistent / Harmonic Balance





