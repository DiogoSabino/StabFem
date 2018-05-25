function [sigma_tab,mode_tab] = SF_FreeMovement_Spectrum(baseflow,sigma_tab,mode_tab,RealShift,ImagShift,STIFFNESS_to_search,mass,nev,stability_analysis)
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
                    mode_tab=[mode_tab em];
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
       
        [ev,em] = SF_Stability(baseflow,'shift',shift,'nev',nev,'type','D','STIFFNESS',STIFFNESS_to_search(1),'MASS',mass,'DAMPING',0,'Frame','R','PlotSpectrum','yes');
        set(gcf, 'Position', get(0, 'Screensize'));%img in fullscream during compute
        disp('TOTO1');
        %[baseflow,em]=SF_Adapt(baseflow,em,'Hmax',10,'InterpError',0.02);
        disp('TOTO2');
        i=1;
        for STIFFNESS=STIFFNESS_to_search
            disp('TOTO3');
            [ev,em] = SF_Stability(baseflow,'shift','cont','nev',nev,'type','D','STIFFNESS',STIFFNESS,'MASS',mass,'DAMPING',0,'Frame','R','PlotSpectrum','yes');
            disp('TOTO4');
            if(mod(i,10)==0)
                disp('TOTO4a');
                [baseflow,em]=SF_Adapt(baseflow,em,'Hmax',10,'InterpError',0.02);
            end
            disp('TOTO5');
            sigma_tab = [sigma_tab ev];
            mode_tab=[mode_tab em];
            i=i+1;
        end
end

end
