function SF_HarmonicForcing(baseflow,Omega_values)
%	File: SF_HarmonicForcing.m
%   Contributours: David Fabre, Diogo Sabino
%   Last Modification: Diogo Sabino, 18 April 2018
global ff ffdir ffdatadir sfdir verbosity ffdataharmonicdir

%% Files management and recover of omega values already calculated
Omega_to_compute=[];

all_data_stored_file=[ffdataharmonicdir, 'Forced_Harmonic2D_Re' num2str(baseflow.Re) 'TOTAL.ff2m'];
if exist(all_data_stored_file)~=0
    %import all data already calculted
    all_data_stored=importFFdata(all_data_stored_file);
    %discover the values of Omega that have not been calculated
    for i=1:size(Omega_values,2)
        if(ismemebertol(Omega_values(i),all_data_stored.OMEGAtab)==0)
            Omega_to_compute=[Omega_to_compute Omega_values(i)]     
        end
    end
end

% Number of omega values to be calculated
NOmega=size(Omega_to_compute,0);

%% .edp command execution
if isempty(Omega_to_compute)==1
    disp('All demanded omegas were already calculated');
    
    
    % And the .edp won't be called
else
    % And the .edp will be called
    % Create the command string:
    stringparam = [Nomega ' '];
    for omega = Omega_to_compute
        stringparam = [stringparam num2str(omega) '  ' ];
    end
    disp(['The following Omega values will be calculated: ' stringparam]);
    
    command = ['echo  ' stringparam ' | ' ff ' FF_Forced_Harmonic_2D_Lateral_Oscillations.edp'];
    error='Error in FF_Forced_Harmonic_2D_Lateral_Oscillations.edp';
    
    mysystem(command,error); %Execute .edp file with FreeFEM++ for specific Re
    
    %Add new omegas calculated
    
    new_data_stored_file=[ffdataharmonicdir 'Forced_Harmonic2D_Re' num2str(baseflow.Re) '.ff2m']);
    new_data_stored=importFFdata(new_data_stored_file);
    
    
    if exist(all_data_stored_file)~=0
        
        for i=1:size(new_data_stored.OMEGAtab,1)
            for j=1:size(all_data_stored.OMEGAtab,1)
                
                if(new_data_stored.OMEGAtab(i)>all_data_stored.OMEGAtab(j))
                    %insert after x3
                    
                end
            end
        end

        %save the data back !!!
        
    else
        % cp all_data_stored_file new_data_stored_file
  
    end
    
    
    
    
end






txt_data = importdata(file);

if isempty(txt_data)==1 %If no value is present
    Omegatab=[];
    Ltab=[];
    Mtab=[]; %Use these infos in DF_HarmonicForcing_Display
    disp('Some error occur. No omega values presented.')
else
    
    %We need to sort the data because .edp file just appends the new omega
    txt_data_ordenated=sortrows(txt_data,1);
    dlmwrite(file,txt_data_ordenated,'delimiter',' ') %refresh the .txt data
    
    %Save data for plot
    Omegatab= txt_data_ordenated(:,1);
    Ltab = txt_data_ordenated(:,3)+1i*txt_data_ordenated(:,4);
    Mtab = txt_data_ordenated(:,5)+1i*txt_data_ordenated(:,6);
    
end

end

%Trash

% % % Unsuccessful essays for don't change stringparam at each iter
% % % numel(num2str(Omega_values))-numel(strfind(num2str(Omega_values),'
% % % '))+numel(Omega_values)*2+4 to know the stringparam length
% % % strlength(stringparam)
