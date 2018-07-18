function [baseflow,sigma_tab] = SF_FreeMovement_Spectrum(stability_analysis,baseflow,sigma_tab,RealShift,ImagShift,STIFFNESS_to_search,mass,nev)
%%
switch stability_analysis
    case 'search' %For search "randomly" in the spectrum
        %RealShift=0.03;
        %ImagShift=0.76;
        %STIFFNESS_to_search=Stiffness_table;
        %sigma_tab=[];
        for STIFFNESS=STIFFNESS_to_search
            for ii=RealShift
                for jj=ImagShift
                    [ev,em] = SF_Stability(baseflow,'shift',ii+jj*1i,'nev',nev,'type','D','STIFFNESS',STIFFNESS,'MASS',mass,'DAMPING',0,'Frame','R','PlotSpectrum','yes');
                    sigma_tab = [sigma_tab ev];
                    %mode_tab=[mode_tab em];It's too heavy in terms of memory
                    set(gcf, 'Position', get(0, 'Screensize'));
                end
            end
        end
        
    case 'modefollow' %For following a mode in the spectrum
        %RealShift=0.03;
        %ImagShift=0.76;
        %STIFFNESS_to_search=Stiffness_table;
        %sigma_tab=[];
        shift=RealShift+ImagShift*1i;
        
        %%%%%%%%%%%plotFF(baseflow,'mesh');
        %Adapt to the first value
        %%%%%%%%%%%[evadapt,emadapt] = SF_Stability(baseflow,'shift',shift,'nev',nev,'type','S','STIFFNESS',STIFFNESS_to_search(1),'MASS',mass,'DAMPING',0,'Frame','R');
        %%%%%%%%%%%[baseflow,em]=SF_Adapt(baseflow,emadapt,'Hmax',1,'InterpError',0.02);
        %%%%%%%%%%%plotFF(baseflow,'mesh');
        
        [ev,em] = SF_Stability(baseflow,'shift',shift,'nev',nev,'type','D','STIFFNESS',STIFFNESS_to_search(1),'MASS',mass,'DAMPING',0,'Frame','R','PlotSpectrum','yes');
        set(gcf, 'Position', get(0, 'Screensize'));%img in fullscream during compute

        
        %for STIFFNESS=STIFFNESS_to_search
        for index_S=1:size(STIFFNESS_to_search,2)
            [ev,em] = SF_Stability(baseflow,'shift','cont','nev',nev,'type','D','STIFFNESS',STIFFNESS_to_search(index_S),'MASS',mass,'DAMPING',0,'Frame','R','PlotSpectrum','yes');
            
            %%%%%%%%%%%if(mod(index_S,10)==0)
            %%%%%%%%%%%    [evadapt,emadapt] = SF_Stability(baseflow,'shift',ev,'nev',nev,'type','S','STIFFNESS',STIFFNESS_to_search(index_S),'MASS',mass,'DAMPING',0,'Frame','R');
            %%%%%%%%%%%    [baseflow,emADAPTED]=SF_Adapt(baseflow,emadapt,'Hmax',1,'InterpError',0.02);
            %%%%%%%%%%%    plotFF(baseflow,'mesh');
            %%%%%%%%%%%end
            sigma_tab = [sigma_tab ev];
            %mode_tab=[mode_tab em];It's too heavy in terms of memory
        end
end

end
