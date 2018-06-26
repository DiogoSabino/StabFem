
function [meanflow,mode] = SF_SelfConsistentDirect(meanflow,mode,varargin)

%%% management of optionnal parameters
    p = inputParser;
   addParameter(p,'Re',meanflow.Re,@isnumeric);
   addParameter(p,'Aguess',-1,@isnumeric);
   addParameter(p,'Fyguess',-1,@isnumeric);
   addParameter(p,'omegaguess',imag(mode.lambda));
   addParameter(p,'sigma',0);
   addParameter(p,'Cyguess',-1,@isnumeric);
   addParameter(p,'Yguess',-1,@isnumeric);
   
   %Parameter for VIV-Diogo
   %addParameter(p,'Yguess',-1,@isnumeric);
   addParameter(p,'STIFFNESS',1000); %its the inverse of U, so it's a Ustar=0
   addParameter(p,'MASS',1000);
   addParameter(p,'DAMPING',0);
   
   parse(p,varargin{:});
   
global ff ffdir ffdatadir sfdir verbosity



if(meanflow.datatype=='BaseFlow')
    disp('### Self Consistent  : with guess from BaseFlow/Eigenmode');
    system(['cp ',ffdatadir, 'BaseFlow.txt ',ffdatadir, 'MeanFlow_guess.txt']);
    system(['cp ',ffdatadir, 'Eigenmode.txt ',ffdatadir, 'SelfConsistentMode_guess.txt']);
    
elseif(meanflow.datatype=='MeanFlow')
    disp('### Self Consistent : with guess from MeanFlow/SCMode');
    system(['cp ',ffdatadir, 'MeanFlow.txt ',ffdatadir, 'MeanFlow_guess.txt']);
    system(['cp ',ffdatadir, 'SelfConsistentMode.txt ',ffdatadir, 'SelfConsistentMode_guess.txt']);

else
    error('wrong type of field for Harmonic balance'); 
end

switch meanflow.mesh.problemtype
    case('2D')
        if(p.Results.Fyguess~=-1)
            disp(['starting with guess Lift force : ' num2str(p.Results.Fyguess) ]);
            solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma)...
                ' L ' num2str(p.Results.Fyguess) ' | ' ff ' '  ffdir 'SelfConsistentDirect_2D.edp'];
        elseif(p.Results.Aguess~=-1)
            disp(['starting with guess amplitude (Energy) ' num2str(p.Results.Aguess) ]);
            solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma)...
                ' E ' num2str(p.Results.Aguess) ' | ' ff ' '  ffdir 'SelfConsistentDirect_2D.edp'];
        else
            solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma)...
                ' none  | ' ff ' '  ffdir 'SelfConsistentDirect_2D.edp'];
            %For comparison with the BIGSPACE
            %disp('Using the BIGSPACE METHOD');
            %solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma)...
            %              ' none  | ' ff ' -v 0 '  ffdir 'SelfConsistentDirect_2D_BIGSPACE.edp'];
        end
    case('2D_VIV')
        disp('DIOGO toto: Entering in VIV HB_O1_VIV.edp');
        
        if(p.Results.Fyguess~=-1)
            disp(['starting with guess Lift force : ' num2str(p.Results.Fyguess) ]);
            solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma) ' L ' num2str(p.Results.Fyguess) ' '...
                num2str(p.Results.STIFFNESS) ' ' num2str(p.Results.MASS) ' ' num2str(p.Results.DAMPING)  ' | ' ff ' '  ffdir 'HB_O1_VIV.edp'];
            
        elseif(p.Results.Aguess~=-1)
            disp(['starting with guess amplitude (Energy) ' num2str(p.Results.Aguess) ]);
            solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma) ' E ' num2str(p.Results.Aguess) ' ' ...
                num2str(p.Results.STIFFNESS) ' ' num2str(p.Results.MASS) ' ' num2str(p.Results.DAMPING)    ' | ' ff ' '  ffdir 'HB_O1_VIV.edp'];
        elseif(p.Results.Yguess~=-1)
            disp(['starting with guess amplitude (Y) ' num2str(p.Results.Aguess) ]);
            solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma) ' Y ' num2str(p.Results.Yguess) ' ' ...
                num2str(p.Results.STIFFNESS) ' ' num2str(p.Results.MASS) ' ' num2str(p.Results.DAMPING)    ' | ' ff ' '  ffdir 'HB_O1_VIV.edp'];

        else
            solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma) ' none ' ...
                num2str(p.Results.STIFFNESS) ' ' num2str(p.Results.MASS) ' ' num2str(p.Results.DAMPING)    ' | ' ff ' '  ffdir 'HB_O1_VIV.edp'];
        end
        
        
        
end

disp('DIOGO TOTO INSIDE DRIVER BEFORE COMPUTATION.')

status = system(solvercommand) %DIOGO: j'avais besoin de verbosity=10 et ce ne marchait pas mysystem...

disp('DIOGO TOTO INSIDE DRIVER AFTER COMPUTATION.')




disp(['#### SELF CONSISTENT CALCULATION COMPLETED with Re = ' num2str(p.Results.Re) ' ; sigma = ' num2str(p.Results.sigma)  ]);
meanflow=importFFdata(meanflow.mesh,'MeanFlow.ff2m');
mode=importFFdata(meanflow.mesh,'SelfConsistentMode.ff2m');

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
