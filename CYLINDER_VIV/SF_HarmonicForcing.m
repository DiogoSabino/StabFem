function [Omegatab,Ltab,Mtab] = SF_HarmonicForcing(baseflow,Omega_values)
%	File: SF_HarmonicForcing.m
%   Contributours: David Fabre, Diogo Sabino
%   Last Modification: Diogo Sabino, 18 April 2018
global ff ffdir ffdatadir sfdir verbosity ffdataharmonicdir

%% Files management and recover of omega values already calculated
file = [ ffdataharmonicdir,'HARMONIC_2D_LateralOscillations_Re', num2str(baseflow.Re),'.txt']
%This next creation isn't need because .edp script append the data and if
% .txt doesn't exist he will created. But for security reasons we put it
if exist(file)==0 %For create the .txt file if doesnt exist;
    mysystem(['touch ' file],'skip');
    disp(['New .txt file HARMONIC_2D_LateralOscillations txt file created for Re=',num2str(baseflow.Re)] )
else
    disp(['.txt file HARMONIC_2D_LateralOscillations already exists for Re=' num2str(baseflow.Re)])
end

old_data_file=importdata(file); %Import omega values already calculated
if isempty(old_data_file)==1 %If no value have been calculated
    Omegas_data_file=[];
else
    Omegas_data_file= old_data_file(:,1);
end

%% Selection only the new omega values
% Next variable isn't needed; created for clarity of the next lines
Omegas_to_finally_comute=Omega_values;

i=1; j=1;
while i<=numel(Omegas_to_finally_comute)
    while j<=numel(Omegas_data_file)
        %Order of next conditions is important, dont change it
        if isempty(Omegas_to_finally_comute)==0 &&i<=numel(Omegas_to_finally_comute) && ismembertol(Omegas_to_finally_comute(i),Omegas_data_file(j))==1
            Omegas_to_finally_comute(i)=[];
        else
            j=j+1;
        end
    end
    i=i+1;
    j=1;
end

Omega_values=Omegas_to_finally_comute; %The useless variable

% From here, Omega_values contains only the values that aren't
% in the .txt file yet, so .edp will append theses ones to the .txt
%% .edp command execution

if isempty(Omega_values)==1
    disp('All demanded omegas were already calculated');
    % And the .edp won't be called
else
    % And the .edp will be called
    % Create the command string:
    stringparam = [];
    for omega = Omega_values;
        stringparam = [stringparam, num2str(omega), '  ' ];
    end
    disp(['The following Omega values will be calculated: ' stringparam]);
    % For stoping .edp internally:
    stringparam = [stringparam, num2str(-1), '  ' ];
    
    command = ['echo  ' ' ', stringparam, ' ' ' | ',ff,' ',' HARMONIC_2D_LateralOscillations.edp'];
    error='Error in HARMONIC_2D_LateralOscillations.edp';
    
    mysystem(command,error); %Execute .edp file with FreeFEM++ for specific Re
end

%% Create the matlab data, to be treated

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
