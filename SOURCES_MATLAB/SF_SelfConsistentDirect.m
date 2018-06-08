
function [meanflow,mode] = SF_SelfConsistentDirect(meanflow,mode,varargin)

%%% management of optionnal parameters
    p = inputParser;
   addParameter(p,'Re',meanflow.Re,@isnumeric);
   addParameter(p,'Aguess',-1,@isnumeric);
   addParameter(p,'Fyguess',-1,@isnumeric); 
   addParameter(p,'omegaguess',imag(mode.lambda));
   addParameter(p,'sigma',0);
   addParameter(p,'Cyguess',-1,@isnumeric);
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
    %DIOGO: j'ai enlevé ces lignes parce que le fichier SF_WNL genere déjà
    %les fichiers guess. Et , en executant SF_WNL avant SF_SC_Direct, si il
    %existe  un meanflow à un Re different, il va efacer le guess qu'on
    %voulez utiliser!!!!!
    %Diogo: neanmoins, ces lignes sont nessessaries pour la partie loop...
    %diogo va revoir son comentaire et repenser sur le sujet. desole .
else
    error('wrong type of field for Harmonic balance'); 
end

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
                  ' None  | ' ff ' '  ffdir 'SelfConsistentDirect_2D.edp'];
    %For comparison with the BIGSPACE
    %disp('Using the BIGSPACE METHOD');
    %solvercommand = ['echo ' num2str(p.Results.Re)  ' ' num2str(p.Results.omegaguess) ' ' num2str(p.Results.sigma)...
    %              ' None  | ' ff ' -v 0 '  ffdir 'SelfConsistentDirect_2D_BIGSPACE.edp'];        
 end
   
   
 status = system(solvercommand) %DIOGO: j'avais besoin de verbosity=10 et ce ne marchait pas mysystem...

               
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
