function Save_Data(Re,m_star,Stiffness_table,U_star,filename,sigma_tab,mode_tab,savedata_dir_version,stability_analysis)
%% 
%Layer Re
Re_dir=['./Final_results_' savedata_dir_version '/Re' num2str(Re) '/'];
if(exist(Re_dir)~=7&&exist(Re_dir)~=5) %je n'est pas compris tres bien cette commande; a voir ensemble apres
    mysystem(['mkdir ' Re_dir]);
end
%Layer mstar
mstar_dir=[Re_dir 'mstar' num2str(m_star) '/'];
if(exist(mstar_dir)~=7&&exist(mstar_dir)~=5) %je n'est pas compris tres bien cette commande; a voir ensemble apres
    mysystem(['mkdir ' mstar_dir]);
end

switch stability_analysis
    case 'spectrum'
        %Save data (change to the good folder manually):
        save([mstar_dir filename '.mat'],'sigma_tab','mode_tab','Stiffness_table','U_star','m_star');
        saveas(gcf,[mstar_dir filename '.fig']);
        saveas(gcf,[mstar_dir filename '.png']);
    case 'grafic'
        saveas(gcf,[mstar_dir filename '.fig']);
        saveas(gcf,[mstar_dir filename '.png']);
    
end

end