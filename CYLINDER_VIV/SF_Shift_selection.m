function [RealShift, ImagShift]=SF_Shift_selection(modename,Re,m_star)


%Manually put the right shift, acordding to the above spectrum(
%For the pair (Re;[m_star])(starting at U_star(1)=3 )
%For the STRUCTURE Mode:
%shift=0+2.1i:
%shift=0+2i:	([60-19];[1000-4.73])
%shift=-0.05+1.7i:                (20;[3])
%shift=-0.1+1.5i:        (20;[2])

%For the FLUID Mode :
%shift=0.05+0.75i:	(60;[20,10,5])
%shift=-0.03+0.75i:               (40,[300-5])
%shift=-0.07+0.62i:                                                                 %%%%(21;[10])
%close all

switch modename{1}
    case('02modeSTRUCTURE')
        if (Re>15&&Re<70)
            
            if(m_star>=0.05&&m_star<0.1)
                   RealShift=-0.07; ImagShift=0.85; %u=1
                   
            elseif(m_star>=0.1&&m_star<0.2)
                RealShift=-0.01; ImagShift=0.68; %u=2

            elseif(m_star>=0.2&&m_star<0.3)
                RealShift=-0.05; ImagShift=0.8; %u=2
                
            elseif(m_star>=0.3&&m_star<0.5)
                RealShift=-0.03; ImagShift=0.76;
                
            elseif(m_star>=0.5&&m_star<1)
                RealShift=-0.05; ImagShift=0.8;
                
            elseif(m_star>=1&&m_star<2&&Re<45)
                RealShift=-0.13; ImagShift=1.1;
            
            elseif(m_star>=1&&m_star<2&&Re>45)
                RealShift=-0.04; ImagShift=1.15;                
                
            elseif(m_star>=2&&m_star<3)
                RealShift=-0.1; ImagShift=1.45;
                
            elseif(m_star>=3&&m_star<4.73)
                RealShift=-0.05; ImagShift=1.7;
                
            elseif(m_star>=4.73&&m_star<20)
                RealShift=0; ImagShift=1.9;
                
            elseif(m_star>=20&&m_star<300)
                RealShift=0; ImagShift=2;
                
            elseif(m_star>=300)
                RealShift=0; ImagShift=2.1;
            end
        else
            disp('Shift not programed');
            RealShift=300; ImagShift=300;
        end
    case('03modeFLUID')
        if (Re>50&&Re<=60)
            RealShift=0.05; ImagShift=0.75;
            
        elseif(Re>45&&Re<=50)
            RealShift=0; ImagShift=0.75;
            
        elseif(Re>40&&Re<=45)
            RealShift=-0.03; ImagShift=0.75;
            
        elseif(Re>35&&Re<=40)
            RealShift=-0.05; ImagShift=0.72;
            
        else
            disp('Shift not programed');
            RealShift=300; ImagShift=300;
        end
         
end

%CHOOSE shift manually:
%RealShift=0; ImagShift=2; %Normally, for Structure
%RealShift=0.04; ImagShift=0.75; %Normally, for Fluid
end