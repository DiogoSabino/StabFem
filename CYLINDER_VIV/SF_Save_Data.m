function SF_Save_Data(save_option,General_data_dir,savedata_dir,Re,m_star,filename,Stiffness_table,U_star,sigma_tab)
%%
%General_data_folder is passed, for the general images

%When 'spectrum' options is used, only one path is allowed and so we need
%to check if it exists and, if not, create it. This justifies the forloop
%nextly. If 'graphic' is used, 2 branchs arose: path is unique or it exists
%multiple path. if path is unique we need the next loop to define it;if
%not, path will store the last this analysed(i think) but nevertheless it
%will not be used.
%So, this justifies why we don't put the next thing in just one of the
%options.

%Checking if directory exists:
%path=[savedata_dir 'Re' num2str(Re) '/mstar' num2str(m_star) '/'];

path={};
%Re_cell=num2cell(Re);
%mstar_cell=num2cell(m_star);

%HERE filename is always just one

for element= 1: size(savedata_dir,2)
    for R=Re
        for m_sta=m_star
            path{end+1}=strcat(savedata_dir{element},'Re',num2str(R),'/mstar',num2str(m_sta),'/');
            if(exist(path{end})~=7&&exist(path{end})~=5) %je n'est pas compris tres bien cette commande; a voir ensemble apres
                %I read in internet that the '-p'(stands for parent) not always work in every shell...
                mysystem(['mkdir -p ' path{end}]);
                disp('Creating new directory for saving data');
            end
        end
    end
end%MUST BE CREATED TO KNOW IF ITS JUST ONLE ON MORE THING THAT ARE PLOTED
%ou entoao passe_se o all path do data treatemetnt e é so fazer o tste

%Choosing what we want to save
switch save_option
    case 'data'
        %Normally, it exists just one cell in path
        %Save data:
        save([path{1} filename{1} '.mat'],'Re','sigma_tab','Stiffness_table','U_star','m_star');
        saveas(gcf,[path{1} filename{1} '.fig']);
        saveas(gcf,[path{1} filename{1} '.png']);
        disp('Data saved!');
    case 'graphic'
        %If we treat just one case:
        if (size(path,2)==1) %which means that only one patth is to be considered
            disp(['Saving image for the selected data with the following name: ' filename{1}]);
            saveas(gcf,[path{1} filename{1} '.fig']);
            saveas(gcf,[path{1} filename{1} '.png']);
        else
            disp('Saving a comparing image of different data');
            %if we compare several cases
            saveas(gcf,[General_data_dir filename{1} '.fig']);
            saveas(gcf,[General_data_dir filename{1} '.png']);

        end
    %case('mesh')%option save mesh
end

end