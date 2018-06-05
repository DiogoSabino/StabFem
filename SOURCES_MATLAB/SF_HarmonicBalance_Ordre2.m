
function [meanflow,mode,mode2] = SF_HarmonicBalance_Ordre2(meanflow,mode,mode2,varargin)

 %%% management of optionnal parameters
     p = inputParser;
    addParameter(p,'Re',meanflow.Re,@isnumeric);
%    addParameter(p,'Aguess',-1,@isnumeric);
%    addParameter(p,'Fyguess',-1,@isnumeric); 
%    addParameter(p,'omegaguess',imag(mode.lambda));
%    addParameter(p,'sigma',0);
%    addParameter(p,'Cyguess',-1,@isnumeric);
    parse(p,varargin{:});


global ff ffdir ffdatadir sfdir verbosity

% if(meanflow.datatype=='BaseFlow')
%     disp('### Self Consistent  : with guess from BaseFlow/Eigenmode');
%     system(['cp ',ffdatadir, 'BaseFlow.txt ',ffdatadir, 'MeanFlow_guess.txt']);
%     system(['cp ',ffdatadir, 'Eigenmode.txt ',ffdatadir, 'SelfConsistentMode_guess.txt']);
    
% % elseif(meanflow.datatype=='MeanFlow')
%    disp('### Self Consistent : with guess from MeanFlow/SCMode');
%    system(['cp ',ffdatadir, 'MeanFlow.txt ',ffdatadir, 'MeanFlow_guess.txt']);
%    system(['cp ',ffdatadir, 'SelfConsistentMode.txt ',ffdatadir, 'SelfConsistentMode_guess.txt']);  
%    system(['cp ',ffdatadir, 'SecondHarmonicMode.txt ',ffdatadir, 'SecondHarmonicMode_guess.txt']);  

% else
%     error('wrong type of field for Harmonic balance'); 
% end

% % if(p.Results.Fyguess~=-1) 
%       disp(['starting with guess Lift force : ' num2str(p.Results.Fyguess) ]);
%      solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma)...
%                   ' L ' num2str(p.Results.Fyguess) ' | ' ff ' '  ffdir 'SelfConsistentDirect_2D.edp'];
%  elseif(p.Results.Aguess~=-1)
%       disp(['starting with guess amplitude (Energy) ' num2str(p.Results.Aguess) ]);
%      solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma)...
%                  ' E ' num2str(p.Results.Aguess) ' | ' ff ' '  ffdir 'SelfConsistentDirect_2D.edp'];
% else
     solvercommand = ['echo ' num2str(p.Results.Re)  ' | ' ff ' '  ffdir 'HarmonicBalance_Order2_2D.edp'];    
% end
   
   
 status = mysystem(solvercommand);

               
disp(['#### HARMOPNIC BALANCE CALCULATION COMPLETED with Re = ' num2str(p.Results.Re) ]);
meanflow=importFFdata(meanflow.mesh,'MeanFlow.ff2m');
mode=importFFdata(meanflow.mesh,'SelfConsistentMode.ff2m');
mode2=importFFdata(meanflow.mesh,'SecondHarmonicMode.ff2m');

if(meanflow.iter<0)
    error('ERROR in SF_HarmonicBalance : Newton iteration did not converge')
end

disp(['#### omega =  ' num2str(imag(mode.lambda)) ]);
%disp(['#### A =  ' num2str(mode.A) ]);


end

%if(nargout>0)
%system('cp MeanFlow.txt MeanFlow_guess.txt');
%system('cp chbase_threshold.txt Self_guess.txt');
%end


%if(nargout>1)
%eigenmode=eigenmodeT;
%system('cp Eigenmode_threshold.txt Eigenmode_guess.txt');
%end
