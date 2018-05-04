function baseflow = SF_Init(meshfile,parameters)
% Matlab/FreeFem driver for generating initial mesh and base flow
%
% usage in one-input mode : baseflow = SF_Init('Mesh.edp')
%
% usage in two-input mode : baseflow = SF_Init('Mesh.edp',params)
%   in this case params is an array containing the parameters needed by the
%   freefem script ;  for instance the dimensions of the mesh  
%
% 'Mesh.edp' must be a FreeFem script which generates a file "mesh.msh",
% a description file "mesh.ff2m", a parameter file "SF_Init.ff2m", and an initial base flow "BaseFlow_init.txt" / "BaseFlow_init.ff2m" 
%
% Version 2.0 by D. Fabre ,  june 2017

global ff ffdir ffdatadir sfdir verbosity

if(exist(ffdatadir)~=7&&exist(ffdatadir)~=5)
    mysystem(['mkdir ' ffdatadir ]); 
else
    mysystem(['rm ' ffdatadir '*.txt ' ffdatadir '*.ff2m ' ffdatadir '*.msh '],'skip');
end

if(exist([ffdatadir 'BASEFLOWS'])~=7)
    mysystem(['mkdir ' ffdatadir 'BASEFLOWS']); 
end
mysystem(['rm ' ffdatadir 'BASEFLOWS/*'],'skip'); 


if(nargin==1)
    command = [ff,' ',meshfile];
else
    stringparam = []; 
    for p = parameters;
        stringparam = [stringparam, num2str(p), '  ' ]; 
    end
    command = ['echo  '' ', stringparam, ' '' | ',ff,' ',meshfile];
end

error = 'ERROR : Freefem not working (path may be wrong, change variable ff in script)';
mysystem(command,error);

   
if(nargout==1)
mesh = importFFmesh('mesh.msh');
mysystem(['cp mesh.msh ' ffdatadir '/mesh_init.msh'],'skip'); 
mysystem(['cp BaseFlow_guess.txt ' ffdatadir 'BASEFLOWS/BaseFlow_init.txt'],'skip'); 
mesh.namefile=[ ffdatadir 'BASEFLOWS/mesh_init.msh'];
baseflow=importFFdata(mesh,'BaseFlow.ff2m');
baseflow.namefile = [ ffdatadir 'BASEFLOWS/BaseFlow_init.txt'];
disp(['      ### INITIAL MESH CREATED WITH np = ',num2str(mesh.np),' points']);

%system(['rm ' ffdatadir 'Eigenmode_guess.txt']);

end
    