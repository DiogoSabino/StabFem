function SF_Data_Treatement(MODE,AXIS,General_data_dir,folder_plot,Re_plot,m_star_plot)
%%
%General_data_folder is passed, with the purpose of being passed to
%SF_Save_Data
%all_paths=[];

all_paths={};
Legend={};
for element= 1: size(folder_plot,2)
    for R=Re_plot
        for m=m_star_plot
            all_paths{end+1}=strcat(folder_plot{element},'Re',num2str(R),'/mstar',num2str(m),'/');
        end
    end
end

%If so pedires pedires dois, mete colres bem. senao, deixao escolher
%color_grafic={'-ob','-or','-og','-oy'};%PUT MORE COLORs
%%%%%%%%%%color_grafic=['-ob';'-or';'-og';'-oy'];  % por a azul sempre em primeiro lugar depois a vermelha...
%tem de ter no minimo o mesmo numero que o maximo que se pode pedir


% % % % % for line= 1: size(folder_plot,1)
% % % % % for Re=Re_plot
% % % % %         for m_star=m_star_plot
% % % % %             all_paths=[all_paths; folder_plot(line,:) 'Re' num2str(Re) '/mstar' num2str(m_star) '/'];
% % % % %         end
% % % % %     end
% % % % % end %For storing all paths

f_subplot1=figure; % Used by all options
%f2=figure; % Used by option 2,3 ...

for element= 1: size(all_paths,2)
    %for line= 1: size(all_paths,1)
    if(exist(all_paths{element})~=7&&exist(all_paths{element})~=5) %je n'est pas compris tres bien cette commande; a voir ensemble apres
        disp('ERROR: The Folder of the data demanded was not calculated yet !!!!');
    else
        switch MODE
            case('Mode:Fluid')
                extracting=[all_paths{element} '03modeFLUID_data.mat'];
                Legend{end+1}=extracting;
                ploting=load( extracting,'Re','m_star','sigma_tab','U_star','Stiffness_table');
                filename={'03modeFLUID'};%For the saving, in the end of this function
            case('Mode:Structure')
                extracting=[all_paths{element} '02modeSTRUCTURE_data.mat'];
                Legend{end+1}=extracting;
                ploting=load( extracting,'Re','m_star','sigma_tab','U_star','Stiffness_table');
                filename={'03modeSTRUCTURE'};%For the saving, in the end of this function
            case('Mode:Both')
                extracting=[all_paths{element} '03modeFLUID_data.mat'];
                Legend{end+1}=extracting;
                extracting2=[all_paths{element} '02modeSTRUCTURE_data.mat'];
                Legend{end+1}=extracting2;
                ploting=load( extracting,'Re','m_star','sigma_tab','U_star','Stiffness_table');
                ploting2=load( extracting2,'Re','m_star','sigma_tab','U_star','Stiffness_table');
                filename={'04modeBOTH'};%For the saving, in the end of this function
        end
        switch AXIS
            case('Axis:sigma_VS_Ustar') %Done
                filename={[filename{1} '_sigma_VS_Ustar']};%For the saving, in the end of this function
                figure(f_subplot1);
                subplot(2,1,1);hold on
                plot(ploting.U_star,real(ploting.sigma_tab),'o-','MarkerSize',2);
                plot(ploting.U_star,ploting.U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
                if(strcmp(MODE,'Mode:Both')==1)
                    plot(ploting2.U_star,real(ploting2.sigma_tab),'o-','MarkerSize',2);
                end
                title('Amplification Rate');
                xlabel('U^*'); ylabel('\lambda_r');
                subplot(2,1,2);hold on;
                plot(ploting.U_star,imag(ploting.sigma_tab),'o-','MarkerSize',2);
                if(strcmp(MODE,'Mode:Both')==1)
                    plot(ploting2.U_star,imag(ploting2.sigma_tab),'o-','MarkerSize',2);
                end
                title('Oscillation Rate');
                xlabel('U^*'); ylabel('\lambda_i');
            case{'Axis:F_LSA_VS_Ustar','Axis:NavroseMittal2016LockInRe60M20','Axis:NavroseMittal2016LockInRe60M5'} %Done %('Axis:F_LSA_VS_Ustar')
                disp('HEREEEEEEEEEE')
                filename={[filename{1} '_F_LSA_VS_Ustar']};%For the saving, in the end of this function
                figure(f_subplot1);
                subplot(2,1,1);hold on
                plot(ploting.U_star,real(ploting.sigma_tab),'o-','MarkerSize',2);
                plot(ploting.U_star,ploting.U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
                if(strcmp(MODE,'Mode:Both')==1)
                    plot(ploting2.U_star,real(ploting2.sigma_tab),'o-','MarkerSize',2);
                end
                title('Amplification Rate');
                xlabel('U^*'); ylabel('\lambda_r');
                subplot(2,1,2);hold on;
                plot(ploting.U_star,imag(ploting.sigma_tab)/(2*pi),'o-','MarkerSize',2);
                if(strcmp(MODE,'Mode:Both')==1)
                    plot(ploting2.U_star,imag(ploting2.sigma_tab)/(2*pi),'o-','MarkerSize',2);
                end
                title('Non-Dimensional Frequency F_{LSA}=\lambda_i/(2\pi)');
                xlabel('U^*'); ylabel('F_{LSA}');
            case{'Axis:fstar_LSA_VS_Ustar','Axis:NavroseMittal2016LockInRe40M10'} %Done
                filename={[filename{1} '_fstar_LSA_VS_Ustar']};%For the saving, in the end of this function
                figure(f_subplot1);
                subplot(2,1,1);hold on
                plot(ploting.U_star,real(ploting.sigma_tab),'o-','MarkerSize',2);
                plot(ploting.U_star,ploting.U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
                if(strcmp(MODE,'Mode:Both')==1)
                    plot(ploting2.U_star,real(ploting2.sigma_tab),'o-','MarkerSize',2);
                end
                title('Amplification Rate');
                xlabel('U^*'); ylabel('\lambda_r');
                subplot(2,1,2);hold on;
                plot(ploting.U_star,imag(ploting.sigma_tab).*ploting.U_star/(2*pi),'o-','MarkerSize',2);
                if(strcmp(MODE,'Mode:Both')==1)
                    plot(ploting2.U_star,imag(ploting2.sigma_tab).*ploting.U_star/(2*pi),'o-','MarkerSize',2);
                end
                title('Frequency Ratio f^*_{LSA}=\lambda_iU^*/(2\pi)');
                xlabel('U^*'); ylabel('f^*_{LSA}');
            case('Axis:sigma_r_VS_Ustar_LSA') %THE SNAKE
                if(strcmp(MODE,'Mode:Both')==1)
                    filename={[filename{1} '_sigma_r_VS_Ustar_LSA']};
                    figure(f_subplot1); hold on;
                    plot((2*pi)./imag(ploting.sigma_tab),real(ploting.sigma_tab),'o-','MarkerSize',2);
                    plot((2*pi)./imag(ploting2.sigma_tab),real(ploting2.sigma_tab),'o-','MarkerSize',2);
                    plot(ploting.U_star,ploting.U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
                    title('Variation of U^*_{LSA}=2\pi/\lambda_r with U^*');
                    xlabel('U^*_{LSA}'); ylabel('\sigma_r');
                else
                    disp('This system of axis must be done with the option: Mode:Both ');
                end
%             case('Axis:NavroseMittal2016LockIn')
%                 filename={[filename{1} '_Comparing_to_Navrose']};
%                 figure(f_subplot1);
%                 subplot(2,1,1);hold on
%                 plot(ploting.U_star,real(ploting.sigma_tab),'o-','MarkerSize',2);
%                 plot(ploting.U_star,ploting.U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
%                 title('Amplification Rate');
%                 xlabel('U^*'); ylabel('\lambda_r');
%                 subplot(2,1,2);hold on;
%                 plot(ploting.U_star,imag(ploting.sigma_tab)/(2*pi),'o-','MarkerSize',2);
%                 title('Non-Dimensional Frequency F_{LSA}=\lambda_i/(2\pi)');
%                 xlabel('U^*'); ylabel('F_{LSA}');
%                 %Column 1: Ustar; %Column 2:Structure MODE; %Column 3: Fluid MODE
%                 Navrose_Re60_mstar20_real = importdata('./Navrose_Data/RE60_M20_real.csv');
%                 Navrose_Re60_mstar20_imag = importdata('./Navrose_Data/RE60_M20_imag.csv');
%                 switch MODE
%                     case('Mode:Structure')
%                         Legend{end+1}='NAVROSE_Structure';
%                         subplot(2,1,1);hold on
%                         plot(Navrose_Re60_mstar20_real.data(:,1),Navrose_Re60_mstar20_real.data(:,3),'--*');
%                         subplot(2,1,2);hold on;
%                         plot(Navrose_Re60_mstar20_imag.data(:,1),Navrose_Re60_mstar20_imag.data(:,3),'--*');
%                     case('Mode:Fluid')
%                                                 Legend{end+1}='NAVROSE_FLUID';
%                         subplot(2,1,1);hold on
%                         plot(Navrose_Re60_mstar20_real.data(:,1),Navrose_Re60_mstar20_real.data(:,2),'--*');
%                         subplot(2,1,2);hold on;
%                         plot(Navrose_Re60_mstar20_imag.data(:,1),Navrose_Re60_mstar20_imag.data(:,2),'--*');
%                     case('Mode:Both')
%                         subplot(2,1,1);hold on
%                         plot(ploting2.U_star,real(ploting2.sigma_tab),'o-','MarkerSize',2);
%                           Legend{end+1}='NAVROSE_stru';
%                             Legend{end+1}='NAVROSE_fluid';
%                         plot(Navrose_Re60_mstar20_real.data(:,1),Navrose_Re60_mstar20_real.data(:,3),'--*');
%                         plot(Navrose_Re60_mstar20_real.data(:,1),Navrose_Re60_mstar20_real.data(:,2),'--*');
%                         subplot(2,1,2);hold on; 
%                         plot(ploting2.U_star,imag(ploting2.sigma_tab)/(2*pi),'o-','MarkerSize',2);
%                         plot(Navrose_Re60_mstar20_imag.data(:,1),Navrose_Re60_mstar20_imag.data(:,3),'--*');
%                         plot(Navrose_Re60_mstar20_imag.data(:,1),Navrose_Re60_mstar20_imag.data(:,2),'--*');
%                        
%                 end
                %...
                %case('spectrum')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'Axis:NavroseMittal2016LockInRe60M20'
%Navrose comparasion, if demanded: 
if(strcmp(AXIS,'Axis:NavroseMittal2016LockInRe60M20')==1||strcmp(AXIS,'Axis:NavroseMittal2016LockInRe60M5')==1||strcmp(AXIS,'Axis:NavroseMittal2016LockInRe40M10')==1)
    switch AXIS
        case('Axis:NavroseMittal2016LockInRe60M20')
            filename={[filename{1} '_Comparing_to_NavroseRe60M20']};
            Navrose_real = importdata('./Navrose_Data/RE60_M20_real.csv');
            Navrose_imag = importdata('./Navrose_Data/RE60_M20_imag.csv');
        case('Axis:NavroseMittal2016LockInRe60M5')
            filename={[filename{1} '_Comparing_to_NavroseRe60M5']};
            Navrose_real = importdata('./Navrose_Data/RE60_M5_real.csv');
            Navrose_imag = importdata('./Navrose_Data/RE60_M5_imag.csv');
        case('Axis:NavroseMittal2016LockInRe40M10')
            filename={[filename{1} '_Comparing_to_NavroseRe40M10']};
            Navrose_real = importdata('./Navrose_Data/RE40_M10_real.csv');
            Navrose_imag = importdata('./Navrose_Data/RE40_M10_imag.csv'); %J AI PAS ENCORE LE FICHIER
    end
    %Column 1: Ustar; %Column 2:Structure MODE; %Column 3: Fluid MODE
    switch MODE
        case('Mode:Structure')
            Legend{end+1}='NavroseMittal Structure';
            subplot(2,1,1);hold on
            plot(Navrose_real.data(:,1),Navrose_real.data(:,3),'--*');
            subplot(2,1,2);hold on;
            plot(Navrose_imag.data(:,1),Navrose_imag.data(:,3),'--*');
        case('Mode:Fluid')
            Legend{end+1}='NavroseMittal Fluid';
            subplot(2,1,1);hold on
            plot(Navrose_real.data(:,1),Navrose_real.data(:,2),'--*');
            subplot(2,1,2);hold on;
            plot(Navrose_imag.data(:,1),Navrose_imag.data(:,2),'--*');
        case('Mode:Both')
            subplot(2,1,1);hold on
            Legend{end+1}='NavroseMittal Structure';
            Legend{end+1}='NavroseMittal Fluid';
            plot(Navrose_real.data(:,1),Navrose_real.data(:,3),'--*');
            plot(Navrose_real.data(:,1),Navrose_real.data(:,2),'--*');
            subplot(2,1,2);hold on;
            plot(Navrose_imag.data(:,1),Navrose_imag.data(:,3),'--*');
            plot(Navrose_imag.data(:,1),Navrose_imag.data(:,2),'--*');   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adding the legend and saving:
%Figure 1
%figure(f_subplot1); hold on
legend(Legend,'Location','southwest','AutoUpdate','off');

if(size(all_paths,2)==1)
    if((strcmp(MODE,'Mode:Fluid')==1||strcmp(MODE,'Mode:Structure')==1)&&strcmp(AXIS,'Axis:sigma_r_VS_Ustar_LSA')==1)
        disp('Image not saved because this system of axis must be done with the option: Mode:Both')
    else
        SF_Save_Data('graphic',General_data_dir,folder_plot,Re_plot,m_star_plot,filename,0,0,0); %Last 3 not used in 'graphic'
        disp('Figure saved !');
    end
else
    disp('Attention: CHOOSE a filename and save the figure, please!!! (run next lines)');
    
end

%Figure 2
%TO DO if needed



end


