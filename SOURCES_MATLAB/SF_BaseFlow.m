function baseflow = SF_BaseFlow(baseflow,varargin) 
% Matlab/SF_ driver for Base flow calculation (Newton iteration)
%
% usage : baseflow = SF_BaseFlow(baseflow1,'Re',Re,[...])
%
% this driver will lanch the "Newton" program of the corresponding
% case. NB if base flow was already created it simply copies it from
% "BASEFLOW" directory.
%
%
% 		NB if for some reason the mesh/baseflow compatibility was lost, use SF__BaseFlow(baseflow,'Re',Re,'type','PREV') 
%	    to recontstruct the structures and reposition the files correctly.
% 
%       similarly to force recomputation even in the case a file exists (for instance just after adaptmesh) use 
%        SF__BaseFlow(baseflow,'Re',Re,'type','NEW')
%
% Version 2.0 by D. Fabre , september 2017
% 

global ff ffdir ffdatadir sfdir verbosity

%%% MANAGEMENT OF PARAMETERS (Re, Mach, Omegax, Porosity...)
% Explanation
% (Mode 1) if parameters are transmitted to the function we use these ones. 
%      (for instance baseflow = SF_BaseFlow(baseflow1,'Re',10)
% (Mode 2) if no parameters are passed and if the field exists in the previous
% baseflow, we take these
%      (for instance SF_BaseFlow(bf) is equivalent to SF_Baseflow(bf,'Re',bf.Re) )
% (Mode 3) if no previous value we will define default values set in the next lines.
%
% This syntax allows to do baseflow=SF_BaseFlow(baseflow) which is useful
% for instance to recompute the baseflow after mesh adaptation.
%
% Parameters currently handled comprise : Re, Omegax, Porosity. 
% usage of Parameter 'type' is to be rationalized...



%% check if fields previously exist (Mode 2) or assign default value (mode 3)
if(isfield(baseflow,'Porosity')) 
    Porosity=baseflow.Porosity; 
else
    Porosity=0; 
end
if(isfield(baseflow,'Omegax')) 
    Omegax=baseflow.Omegax; 
else
    Omegax = 0; 
end

%%% check which parameters are transmitted to varargin (Mode 1) 
    p = inputParser;
   addParameter(p,'Re',baseflow.Re,@isnumeric); % Reynolds
   addParameter(p,'Omegax',Omegax,@isnumeric); % rotation rate (for swirling body)
   addParameter(p,'Porosity',Porosity,@isnumeric); % For porous body
   addParameter(p,'type','Normal',@ischar); % mode 
   parse(p,varargin{:});
   
% Now the right parameters are in p.Results   
   Re = p.Results.Re;
   Omegax = p.Results.Omegax;
   Porosity=p.Results.Porosity;


%% SELECTION OF THE SOLVER TO BE USED DEPENDING ON THE CASE

switch(baseflow.mesh.problemtype)
    
    case ('AxiXR')
        % Newton calculation for axisymmetric base flow
        if(verbosity>1)  disp('## solving base flow (axisymmetric case)'); end
        solvercommand = ['echo ' num2str(Re) ' | ',ff,' ',ffdir,'Newton_Axi.edp'];
        
    case('AxiXRPOROUS') % axisymmetric WITH SWIRL
        if(verbosity>1)  disp('## solving base flow (axisymmetric case WITH SWIRL)'); end
        solvercommand = ['echo ' num2str(Re) ' ' num2str(p.Results.Omegax) ' ' num2str(p.Results.Porosity) ' | ',ff,' ',ffdir,'Newton_AxiSWIRL.edp']
        
    case('2D')
        if(verbosity>1)  disp('## solving base flow (2D CASE)'); end
        solvercommand = ['echo ' num2str(Re) ' | ',ff,' ',ffdir,'Newton_2D.edp'];
    case('2D_VIV')
        if(verbosity>1)  disp('## solving base flow (2D_VIV CASE=2D)'); end
        solvercommand = ['echo ' num2str(Re) ' | ',ff,' ',ffdir,'Newton_2D.edp'];
        
   % case (other cases...)
        
end %switch

error = 'ERROR : SF_ base flow computation aborted';
        
%% Selection of what to do according to the parameters
if (strcmp(p.Results.type,'PREV')==1)
	% recover base flow from previous adapted case 
    disp(['      ### FUNCTION SF_BaseFlow : recovering previous adapted mesh/baseflow for Re = ', num2str(Re)]);
	file = [ ffdatadir '/BASEFLOWS/BaseFlow_adapt_Re' num2str(Re) '.txt ' ];
         system(['cp ' file ffdatadir ' BaseFlow_guess.txt']);
    file = [ ffdatadir '/BASEFLOWS/mesh_adapt_Re' num2str(Re) '.msh ' ];
         system(['cp ' file  ffdatadir ' mesh.msh']);
    mysystem(solvercommand,error); %needed to generate .ff2m file
    mesh = importFFmesh('mesh.msh');
    mesh.namefile=[ffdatadir '/BASEFLOWS/mesh_adapt_Re' num2str(baseflow.Re) '.msh'];
    baseflow.mesh=mesh;
    baseflow = importFFdata(baseflow.mesh,'BaseFlow.ff2m'); 
    baseflow.namefile = [ ffdatadir 'BASEFLOWS/BaseFlow_Re' num2str(Re) '.txt'];
    baseflow.iter=0;
    
elseif(exist([ ffdatadir '/BASEFLOWS/BaseFlow_Re' num2str(Re) '.txt'])==2&&strcmp(p.Results.type,'NEW')~=1)
    %
	disp(['FUNCTION SF_BaseFlow : base flow already computed for Re = ', num2str(Re)]);
	system(['cp ' ffdatadir '/BASEFLOWS/BaseFlow_Re' num2str(Re) '.txt  ' ffdatadir 'BaseFlow.txt']);
	system(['cp ' ffdatadir '/BASEFLOWS/BaseFlow_Re' num2str(Re) '.txt  ' ffdatadir 'BaseFlow_guess.txt']);
	system(['cp ' ffdatadir '/BASEFLOWS/BaseFlow_Re' num2str(Re) '.ff2m ' ffdatadir 'BaseFlow.ff2m']);
	baseflow = importFFdata(baseflow.mesh,[ffdatadir 'BaseFlow.ff2m']); 
    baseflow.namefile = [ ffdatadir 'BASEFLOWS/BaseFlow_Re' num2str(Re) '.txt'];
	baseflow.iter=0;
        
else
%% Do the baseflow calculation
 
    if(verbosity>0)
        disp(['      ### FUNCTION SF_BaseFlow : computing base flow for Re = ', num2str(Re)]);
    end
    %system(['cp ' baseflow.namefile ' BaseFlow_guess.txt']);

    mysystem(solvercommand,error);  %%% CALL NEWTON SOLVER

        
    if(exist([ffdatadir,'BaseFlow.txt'])==0)%Because if was deleted by Newton.edp
        error('ERROR : SF_ base flow computation did not converge');
    end
        
	system(['cp ' ffdatadir 'BaseFlow.txt ' ffdatadir 'BASEFLOWS/BaseFlow_Re' num2str(Re) '.txt']);
	system(['cp ' ffdatadir 'BaseFlow.ff2m ' ffdatadir 'BASEFLOWS/BaseFlow_Re' num2str(Re) '.ff2m']);
	baseflow = importFFdata(baseflow.mesh,'BaseFlow.ff2m'); 
    baseflow.namefile = [ ffdatadir 'BASEFLOWS/BaseFlow_Re' num2str(Re) '.txt'];
        
    %system(['cp BaseFlow.txt BaseFlow_guess.txt']);    
 end

 
              
 
 
if(baseflow.iter>=1)
    message = ['      # Base flow converged in ',num2str(baseflow.iter),' iterations '];
    if(isfield(baseflow,'Drag')==1) %% adding drag information for blunt-body wake
        message = [message , '; Drag = ',num2str(baseflow.Drag)];
    end    
    if(isfield(baseflow,'Lx')==1) %% adding drag information for blunt-body wake
        message = [message , '; Lx = ',num2str(baseflow.Lx)];
    end
    if(isfield(baseflow,'deltaP0')==1) %% adding pressure drop information for jet flow
        message = [message , '; deltaP0 = ',num2str(baseflow.deltaP0)];
    end
 disp(message);
else
 disp(['      ### Base flow recovered from previous computation for Re = ' num2str(Re)]);   
end
%if(nargout==1)
%baseflow = importFFdata(baseflow.mesh,['BaseFlow.ff2m']);
%baseflow.namefile = [ ffdatadir '/BaseFlow/BaseFlow_Re' num2str(Re) '.txt'];
%end

end

