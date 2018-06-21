
function [baseflow,eigenmode] = SF_FindThreshold(baseflow,eigenmode)

% Direct computation of instability threshold

global ff ffdir ffdatadir sfdir verbosity

system(['cp ' ffdatadir 'Eigenmode.txt ' ffdatadir 'Eigenmode_guess.txt']);

switch baseflow.mesh.problemtype
    case('2D') %For CILINDER
        disp('Executing FindThreshold2D.edp.');
        solvercommand = [ff ' ' ffdir 'FindThreshold2D.edp'];
        status = mysystem(solvercommand);
        disp('Execution of FindThreshold2D.edp done.');
    case('2D_VIV') %For CILINDER_VIV
        disp('Executing FindThreshold2D_VIV.edp.');
        solvercommand = [ff ' ' ffdir 'FindThreshold2D_VIV.edp'];
        status = mysystem(solvercommand);
        disp('Execution of FindThreshold2D_VIV.edp done.');  
end

        
disp(['#### Direct computation of instability threshold ']);

baseflowT=importFFdata(baseflow.mesh,'BaseFlow_threshold.ff2m');
eigenmodeT=importFFdata(baseflow.mesh,'Eigenmode_threshold.ff2m');


disp(['#### Re_c =  ' num2str(baseflowT.Re) ]);
disp(['#### lambda_c =  ' num2str(eigenmodeT.lambda) ]);

if(nargout>0)
baseflow=baseflowT;
system(['cp ',ffdatadir,'BaseFlow_threshold.txt ',ffdatadir,'BaseFlow.txt']);
system(['cp ',ffdatadir,'BaseFlow_threshold.txt ',ffdatadir,'BaseFlow_guess.txt']);
system(['cp ',ffdatadir,'BaseFlow_threshold.ff2m ',ffdatadir,'BaseFlow.ff2m']);
system(['cp ',ffdatadir,'BaseFlow_threshold.ff2m ',ffdatadir,'BaseFlow_guess.ff2m']);
end


if(nargout>1)
eigenmode=eigenmodeT;
system(['cp ',ffdatadir,'Eigenmode_threshold.txt ',ffdatadir,'Eigenmode_guess.txt']);
system(['cp ',ffdatadir,'Eigenmode_threshold.txt ',ffdatadir,'Eigenmode.txt']);
system(['cp ',ffdatadir,'Eigenmode_threshold.ff2m ',ffdatadir,'Eigenmode_guess.ff2m']);
system(['cp ',ffdatadir,'Eigenmode_threshold.ff2m ',ffdatadir,'Eigenmode.ff2m']);
end



        