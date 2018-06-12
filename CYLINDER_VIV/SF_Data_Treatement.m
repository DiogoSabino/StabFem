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
            case('Axis:sigma_VS_Ustar')
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
            case('Axis:F_LSA_VS_Ustar')
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
            case('Axis:f_LSA_VS_Ustar')
                filename={[filename{1} '_f_LSA_VS_Ustar']};%For the saving, in the end of this function
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
                title('Frequency Ratio f_{LSA}=\lambda_iU^*/(2\pi)');
                xlabel('U^*'); ylabel('f^*_{LSA}');
                
                
        end
        % if(strcmp(AXIS,'sigma_r_VS_Ustar_LSA')==1)
        %continuar em casa
        % % % % % %         switch case_to_treat %This function could be better design, but...
        % % % % % %             case('Mode:Fluid;Axis:sigma_VS_Ustar')
        % % % % % %                 %Real and imaginary sigma values in function of Ustar
        % % % % % %                 extracting=[all_paths{element} '03modeFLUID_spectrum.mat']; %it's a string because I didnt put {}
        % % % % % %                 ploting=load( extracting,'Re','m_star','sigma_tab','U_star','Stiffness_table');
        % % % % % %                 figure(f_subplot1);
        % % % % % %                 subplot(2,1,1);hold on
        % % % % % %                 plot(ploting.U_star,real(ploting.sigma_tab),color_grafic{element},'MarkerSize',2);
        % % % % % %                 plot(ploting.U_star,ploting.U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
        % % % % % %                 title('Amplification Rate');
        % % % % % %                 xlabel('U^*'); ylabel('\sigma_r');
        % % % % % %                 subplot(2,1,2);hold on;
        % % % % % %                 plot(ploting.U_star,imag(ploting.sigma_tab),color_grafic{element},'MarkerSize',2);
        % % % % % %                 title('Oscillation Rate');
        % % % % % %                 xlabel('U^*'); ylabel('\sigma_i');
        % % % % % %                 filename={'03modeFLUI'};%For the saving, in the end of this function
        % % % % % %             case('Mode:Structure;Axis:sigma_VS_Ustar')
        
        
        %case('Mode:Fluid;Navrose')
        %para importar data do navrose
        %Navrose_Re60_mstar20 = importdata('./Navrose_Data/RE60_M20_real.csv'); %meter os bons nomes
        %...
        %case('spectrum')
    end
    
    
end
%i=i+1; %For iterating in color

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adding the legend and saving:
%Figure 1 refazer por causa do modo
%figure(f_subplot1); hold on
legend(Legend,'Location','southwest','AutoUpdate','off');

if(size(all_paths,2)==1)
    SF_Save_Data('graphic',General_data_dir,folder_plot,Re_plot,m_star_plot,filename,0,0,0); %Last 3 not used in 'graphic'
    disp('Figure saved !');
else
    disp('Attention: CHOOSE a filename and save the figure, please!!! (run next lines)');
    
end

%Figure 2
%TO DO if needed



end


