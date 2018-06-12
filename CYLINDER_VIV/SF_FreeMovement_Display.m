function SF_FreeMovement_Display(U_star,sigma_tab,which_plot,modename)
%%
figure;hold on;
set(gcf, 'Position', get(0, 'Screensize'));
color_grafic=['-or';'-ob';'-og';'-oy']; %PUT MORE COLORs

switch which_plot
    case 'sigma_VS_Ustar' %FEITO
        subplot(2,1,1);hold on
        plot(U_star,real(sigma_tab),'-or','MarkerSize',2);
        plot(U_star,U_star*0,'--k','LineWidth',0.1)
        title('Amplification Rate');
        legend(modename)
        xlabel('U^*'); ylabel('\sigma_r');
        
        subplot(2,1,2);hold on;
        plot(U_star,imag(sigma_tab),'-or','MarkerSize',2);
        title('Oscillation Rate');
        legend(modename)
        xlabel('U^*'); ylabel('\sigma_i');
        
        % save data outside this function
    case 'F_LSA_VS_Ustar' %FEITO
        i=1;
        while(i<=size(sigma_tab,1))
            subplot(2,1,1);hold on
            plot(U_star,real(sigma_tab(i,:)),color_grafic(i,:),'MarkerSize',2);
            plot(U_star,U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
            title('Amplification Rate');
            xlabel('U^*'); ylabel('\sigma_r');
            
            subplot(2,1,2);hold on;
            plot(U_star,imag(sigma_tab(i,:))/(2*pi),color_grafic(i,:),'MarkerSize',2);
            title('Non-Dimensional Frequency F_{LSA}=\sigma_i/(2\pi)');
            xlabel('U^*'); ylabel('F_{LSA}');
            
            i=i+1;
        end
        subplot(2,1,1);legend(modename); subplot(2,1,2);legend(modename);
        % save data outside this function
    case 'f_LSA_VS_Ustar' %FEITO
        i=1;
        while(i<=size(sigma_tab,1))
            subplot(2,1,1);hold on
            plot(U_star,real(sigma_tab(i,:)),color_grafic(i,:),'MarkerSize',2);
            plot(U_star,U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
            title('Amplification Rate');
            %legend(modename(i,:))
            xlabel('U^*'); ylabel('\sigma_r');
            
            subplot(2,1,2);hold on;
            plot(U_star,imag(sigma_tab(i,:))/2/pi.*U_star,color_grafic(i,:),'MarkerSize',2);
            title('Oscillation Rate');
            %legend(modename(i,:))
            xlabel('U^*'); ylabel('f^*_{LSA}');
            i=i+1;
        end
        subplot(2,1,1);legend(modename); subplot(2,1,2);legend(modename);
        % save data outside this function
    case 'Ustar_LSA_VS_sigma_i'
        i=1;
        while(i<=size(sigma_tab,1))
            plot((2*pi)./imag(sigma_tab(i,:)),real(sigma_tab(i,:)),color_grafic(i,:),'MarkerSize',2);
            plot(U_star,U_star*0,'--k','LineWidth',0.1,'HandleVisibility','off')
            title('Image from p.574, Navrose');
            %legend(modename(i,:))
            xlabel('U^*_{LSA}'); ylabel('\sigma_r');
            
            i=i+1;
        end
        legend(modename);
end

end