%***********************************************************%
%                 >> IN THE NAME OF ALLAH <<                %
%Moment-Curvature Analysis of Steel Box and Concrete section%
%-----------------------------------------------------------%
%    This program is written by Salar Delavar Ghashghaei    %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given:Section Properties , Concrete properties ,           %
% Reinforcing steel properties                              %
%Calculate: Moment-Curavture                                %
% Note: No limit for accounting plurality steel rebar       %
% Newton-Raphson Method : Tangent procedure                 %
%***********************************************************%
%      ____________________________________________         %
%   _ |   ______________________________________   |        %
%   | |  |                                      |  |        %
%     |  |     #     #     #     #    #    #    |  |        %
%     |  |     #                           #    |  |        %
%   b |  |    As1   As2   As3   As4  As5  As6   |  |        %
%     |  |     #                           #    |  |        %
%   | |  |     #     #     #     #    #    #    |  |        %
%   _ |  |______________________________________|  |        %
%     |____________________________________________|        %
%        |<-                 h                ->|tp|        %
%       |<-d1->|                                            %
%       |<-  d2   ->|                                       %
%       |<-     d3      ->|                                 %
%       |<-        d4          ->|                          %
%       |<-            d5          ->|                      %
%       |<-               d6             >|                 %
%    X                                                      %
%    ^                                                      %
%    |             (Moment - Curvature along X axis)        %
%    |                                                      %
%    +----> Y                                               %
%***********************************************************%
clear all;close all;clc;
%% Section Properties
b=500;% [mm]
h=500;% [mm]
%As:  As1      As2
As=[2454.29 2454.29]; % NOTE: As1 & As2 = 5fi25
%d:d1  d2 
d=[50 450];
tp=10;% [mm] Box section thickness on Top and Bottom Plate
hb=2*tp+h;% [mm] Height of Section
wb=2*tp+b;% [mm] Width of Section
%% Steel Section Properties
fy =240;% [N/mm^2] Yield strength of steel section
Es =2e5;% [N/mm^2] Modulus of elasticity of steel section
fu=1.5*fy;% Ultimate steel stress
ey=fy/Es;% Yield steel strain
esh=0.025;% Strain at steel strain-hardening
esu=0.35;% Ultimate steel strain
Esh=(fu-fy)/(esu-esh);
%% Concrete Properties
fc =25;% [N/mm^2] Unconfined concrete strength
ecu=0.004;% Ultimate concrete strain
Ec=5000*sqrt(fc);
ec0=(2*fc)/Ec;
fct=-0.7*sqrt(fc);% Concrete tension stress
ect1=(2*fct)/Ec;ect2=(2.625*fct)/Ec;ect3=(9.292*fct)/Ec;% Concrete tension strain
%% Reinforcing steel Properties
fys =400;% Yield strength of reinforcing steel (N/mm^2)
Ess =2e5;% Modulus of elasticity of steel (N/mm^2)
fus=1.5*fys;% Ultimate steel stress
eys=fys/Ess;% Yeild steel strain
eshs=0.01;% Strain at steel strain-hardening
esus=0.09;% Ultimate steel strain
Eshs=(fus-fys)/(esus-eshs);
N=1000;% Number of steel section Fiber 
itermax = 4000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
x=.25*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure ------------------------%
N1=1*N/10;% Number of steel section top flange Fiber  
N2=8*N/10;% Number of steel section top flange Fiber  
N3=1*N/10;% Number of steel section  web Fiber 
An=size(As,2);R1=(1/N1);R2=(1/N2);R3=(1/N3);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*tp;end % distance of each top steel plate Fiber
for k=1:N2;c2(k)=tp+(.5*R2+R2*(k-1))*h;end % distance of each steel and concrete section Fiber 
for k=1:N3;c3(k)=tp+h+(.5*R3+R3*(k-1))*tp;end % distance of each bottom steel plate Fiber
c=[c1 c2 c3];
ES=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q=size(ES,2);
% ES=ey:.01:esu;q=size(ES,2);
for j=1:q;
    eS=ES(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        %---------------- As -------------%
        for u=1:An;
        ess=eS*(x-(tp+d(u)))/x;
       if and(ess>=0,ess<=eys)
        fsa(u)=Ess*ess;
        fstana(u)=(Ess*eS*(tp+d(u)))/(x)^2;
       elseif and(ess<0,ess>=(-eys))
        fsa(u)=Ess*ess;
        fstana(u)=(Ess*eS*(tp+d(u)))/(x)^2;
       elseif  and(ess>eys,ess<eshs)
        fsa(u)=fys;
        fstana(u)=0;
       elseif and(ess<=(-eys),ess>(-eshs))
        fsa(u)=-fys;
        fstana(u)=0;
       elseif  and(ess>=eshs,ess<esus)
        fsa(u)=fus-(fus-fys)*((esus-abs(ess))/(esus-eshs))^2;
        fstana(u)=(2*eS*(tp+d(u))*(fus-fys)*(((eS*(tp+d(u)))/x)+esus-eS))/(x^2*(esus-eshs)^2);
       elseif and(ess<=(-eshs),ess>(-esus))
        fsa(u)=-fus+(fus-fys)*((esus-abs(ess))/(esus-eshs))^2;
        fstana(u)=(2*eS*(tp+d(u))*(fus-fys)*(((eS*(tp+d(u)))/x)+esus-eS))/(x^2*(esus-eshs)^2);
       elseif or(ess>=esus,ess<=(-esus))
        fsa(u)=0;
        fstana(u)=0;
       end
        FSA(j,u)=fsa(u);
        FsA(u)=As(u)*fsa(u);
        FstanA(u)=As(u)*fstana(u);% tangent steel force
        end
        
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c(z))/x;       
        if and(z>=1,z<=roundn(1*N/10,0)) % TOP PLATE STEEL SECTION    
        %---------------- Fs -------------%
       if and(es>=0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=wb*R1*tp*fs;% force on top flange
     Fstan(z)=wb*R1*tp*fstan;% tangent steel force
        elseif and(z>=roundn(1*N/10,0),z<=roundn(9*N/10,0)) % TOP CONCRETE SECTION AND WEB PLATE STEEL SECTION
            
        if and(es>=0,es<ec0)
        fsc=fc*((2*es/ec0)-(es/ec0)^2);
        fstanc=((2*fc)/(ec0^2*x^3))*(ec0*eS*c(z)*x-c(z)*eS^2*x+2*eS^2*c(z)^2);
      elseif and(es>=ec0,es<ecu)
        fsc=fc*(1-(0.15*(es-ec0)/(ecu-ec0)));
        fstanc=-(3*eS*c(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif es >= ecu
         fsc=0; 
         fstanc=0;
      elseif and(es<0,es>=ect1)
        fsc=0.5*Ec*es;
        fstanc=(0.5*Ec*c(z)*eS)/x^2;
      elseif and(es<ect1,es>=ect2)
        fsc=fct-(0.5*fct/(ect2-ect1))*(es-ect1);
        fstanc=-(0.5*fct*eS*c(z))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        fsc=.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2);
        fstanc=-(0.5*fct*eS*c(z))/((ect3-ect2)*x^2);
      elseif es<ect3
        fsc=0;
        fstanc=0;
        end 
      
       if and(es>=0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=b*R2*h*fsc+tp*2*R2*h*fs;% force on web
     Fstan(z)=b*R2*h*fstanc+tp*2*R2*h*fstan;% tangent steel force
        elseif and(z>roundn(9*N/10,0),z<=N) % BOTTOM PLATE STEEL SECTION 
       if and(es>=0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=wb*R3*tp*fs;% force on bottom flange
      Fstan(z)=wb*R3*tp*fstan;% tangent steel force 
        end 
     Ss(z)=es;SS(j,z)=Ss(z);% sction Fiber Strain
     CFS(j,z)=fs;% sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);FsATOTAL=sum(FsA);A=FsTOTAL+FsATOTAL;
    FsTOTAL_tan=sum(Fstan);FsATOTAL_tan=sum(FstanA);A_tan=FsTOTAL_tan+FsATOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eS,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    for u=1:An;e(u)=x-(tp+d(u));end % distance of each rebar from neutral axis        
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc(k)=x-c(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc'+FsA*e')*10^-6)
        end                
TSCSs(j)=eS*(x-(tp+d(1)))/x;TSCS(j)=fsa(1);% Top Steel compression strain-stress
BSCSs(j)=eS*(x-(tp+d(2)))/x;BSCS(j)=fsa(2);% Bottom Steel compression strain-stress       
eST=eS*(x-tp)/x;if and(eST>0,eST<ec0);ECC(j)=eST;TUCS(j)=fc*((2*eST/ec0)-(eST/ec0)^2);elseif and(eST>=ec0,eST<ecu);ECC(j)=eST;TUCS(j)=fc*(1-(0.15*(eST-ec0)/(ecu-ec0)));end% Top Unconfined compression stress         
if and(eS>0,eS<ey);EST(j)=eS;TUST(j)=Es*eS;elseif and(eS>=ey,eS<esh);EST(j)=eS;TUST(j)=fy;elseif and(eS>=esh,eS<esu);EST(j)=eS;TUST(j)=fu-(fu-fy)*((esu-abs(eS))/(esu-esh))^2;elseif and(eS<0,eS>(-ey));EST(j)=eS;TUST(j)=Es*eS;elseif and(eS<=(-ey),eS>(-esh));EST(j)=eS;TUST(j)=-fy;elseif and(eS<=(-esh),eS>(-esu));EST(j)=eS;TUST(j)=-fu+(fu-fy)*((esu-abs(eS))/(esu-esh))^2;end% Top Steel section compression stress
eSB=eS*(x-hb)/x;if and(eSB>0,eSB<ey);ESB(j)=eSB;TUSB(j)=Es*eSB;elseif and(eSB>=ey,eSB<esh);ESB(j)=eSB;TUSB(j)=fy;elseif and(eSB>=esh,eSB<esu);ESB(j)=eSB;TUSB(j)=fu-(fu-fy)*((esu-abs(eSB))/(esu-esh))^2;elseif and(eSB<0,eSB>(-ey));ESB(j)=eSB;TUSB(j)=Es*eSB;elseif and(eSB<=(-ey),eSB>(-esh));ESB(j)=eSB;TUSB(j)=-fy;elseif and(eSB<=(-esh),eSB>(-esu));ESB(j)=eSB;TUSB(j)=-fu+(fu-fy)*((esu-abs(eSB))/(esu-esh))^2;end% Bottom Steel section compression stress
% Calculate Moment and Curavture
Cur(j)=(eS/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*cc'+FsA*e')*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
Cur=[0 Cur];Mom=[0 Mom];
s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));end % Flextural Rigidity
if eS == ecu;fprintf('\n      ##  Steel Plate Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
%% Steel Section bilinear fitting
SIZE=size(Mom,2);
for i=1:SIZE-1;
    hh(i) = Cur(i+1)-Cur(i);
    Aa(i)=(Mom(i)+Mom(i+1))*0.5*hh(i);
end
Area=sum(Aa);k0 =Mom(2)/Cur(2);
fiy = (Mom(i+1)*max(Cur)*0.5-Area)/(Mom(i+1)*0.5 - k0*max(Cur)*0.5);
My = k0*fiy;
X = [0 fiy max(Cur)];Y = [0 My Mom(i+1)];
disp('+=============================+');
disp('= Steel Section curve fitted  =');
disp('     Curvature    Moment       ');
disp('       (1/m)      (kN.m)       ');
disp('-------------------------------');
disp([X' Y']);
disp('+==============================+');
%% EI and Ductility_Rito of Steel Section
Elastic_EI=Y(2)/X(2);
Plastic_EI=(Y(3)-Y(2))/(X(3)-X(2));
Ductility_Rito=X(3)/X(2);
Over_Strength_Factor=Y(3)/Y(2);
Material_Ductility_Rito=esu/esh;
fprintf('+------------------------------------------+\n')
fprintf(' Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI)
fprintf(' Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI)
fprintf(' Steel Material Ductility Ratio : %5.2f\n',Material_Ductility_Rito)
fprintf(' Steel Section Ductility Ratio : %5.2f\n',Ductility_Rito)
fprintf(' Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor)
fprintf('+------------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% DATA
Cur01=[         0
    0.0010
    0.0021
    0.0032
    0.0044
    0.0060
    0.0347
    0.0530
    0.0710
    0.0898
    0.1080
    0.2791
    0.5489
    0.8173
    1.0868
       NaN];
Mom01=10^3*[         0
    0.3062
    0.5841
    0.8598
    1.1134
    1.2863
    1.4702
    1.3990
    1.3815
    1.3876
    1.4015
    1.5953
    1.1762
    1.2765
    1.3426
       NaN];
%% SAP2000 - Moment-Curvature of C50x50-10fi25-T10
CurP=[0
    5.425*10^-3
    6.961*10^-3
    1.249*10^-2
    4.045*10^-2
    5.09*10^-2
    5.203*10^-2
    6.682*10^-2
    9.197*10^-2
    4.47*10^-1
    4.648*10^-1
    6.721*10^-1
    8.328*10^-1
    9.692*10^-1
    1.19
    1.35];
MomP=[0
    1219.02
    1328.73
    1441.96
    1466.03
    1357.04
    1347.87
    1332.36
    1330.3
    1711.522
    1136.61
    1220.8
    1269.67
    1298.4
    1328.47
    1340.93];
%% imaging
figure(1)
IMAGE1=imread('FiberCompositeBoxSteelSectionMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE1=imread('FiberCompositeBoxSteelSectionMomentCurvature-image2.jpg');
image(IMAGE1);axis image;axis off;
figure(3)
IMAGE2=imread('FiberCompositeSectionMomentCurvature-image1.jpg');
image(IMAGE2);axis image;axis off;
%% Plot
figure(4)
plot(TSCSs,TSCS,BSCSs,BSCS,'g--','LineWidth',3);xlabel('Steel strain (\epsilon_s)');ylabel('Steel stress (N/mm^2)')
title(' Plotting of the Top  and Bottom steel rebar strain-stress during the analysis','Color','b');
legend('Top Steel Rebar','Bottom Steel Rebar','Location','NorthEastOutside');grid on;
figure(5)
P1=plot(EST,TUST,ESB,TUSB,'r--');set(P1,'LineWidth',3);
xlabel('Strain');ylabel('Stress (kN/mm^2)')
legend('Top Steel Fiber Section','Bottom Steel Fiber Section','Location','NorthEastOutside');grid on;
title('Top and Bottom Steel Fiber Section Stress-Strain Diagram','FontSize',12,'color','b')
figure(6)
P1=plot(ECC,TUCS);set(P1,'LineWidth',3);
xlabel('Strain');ylabel('Stress (kN/mm^2)');grid on;
title('Top Concrete Fiber Section Stress-Strain Diagram','FontSize',12,'color','b')
figure(7)
p1=plot(I,DX,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
figure(8)
p1=plot(I,IT,'b--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
figure(9)
plot(ES,XX,'--','LineWidth',3);xlabel('Top Concrete strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top Steel strain and Neuteral axis ','Color','b');grid on
figure(10)
plot(ES,CUR,'--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top Steel strain and Curvature ','Color','b');grid on
figure(11)
plot(SS(1,:),(hb-c),SS(roundn(.25*q,0),:),(hb-c),SS(roundn(.5*q,0),:),(hb-c),SS(roundn(.75*q,0),:),(hb-c),SS(roundn(q,0),:),(hb-c),'g--','LineWidth',3);xlabel('Strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the steel and concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(12)
plot(CFS(1,:),(hb-c),CFS(roundn(.25*q,0),:),(hb-c),CFS(roundn(.5*q,0),:),(hb-c),CFS(roundn(.75*q,0),:),(hb-c),CFS(roundn(q,0),:),(hb-c),'g--','LineWidth',3);xlabel('Stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the steel and concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(13)
plot(CUR,EI,'black--','LineWidth',3)
title('#  Flextural Rigidity(EI)-CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(14)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral axis-Flextural Rigidity(EI) Diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(15)
plot(Cur,Mom,X,Y,'r--','LineWidth',3)
title(['# COMPOSITE SECTION MOMENT-CURVATURE DIAGRAM #','  EI_e_l_a : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Composite','Composite - bilinear fitted','Location','NorthEastOutside');grid on;
figure(16)
plot(Cur01,Mom01,CurP,MomP,'r--','LineWidth',3)
title(['# COMPOSITE SECTION SECTION MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Analysis','SAP2000','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberCompositeBoxSteelSectionMomentCurvature-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*Moment-Curvature Analysis of Steel Box and Concrete section*\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*     This program is written by salar delavar ghashghaei   *\n');  
fprintf(fid,'*          E-mail:salar.d.ghashghaei@gmail.com              *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*Unit: Newton-Milimeter                                     *\n');
fprintf(fid,'*Given:Section Properties , Concrete properties ,           *\n');
fprintf(fid,'* Reinforcing steel properties                              *\n');
fprintf(fid,'*Calculate: Moment-Curavture                                *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*   _    ______________________________________             *\n');
fprintf(fid,'*   |   |                                      |            *\n');
fprintf(fid,'*       |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*   b   |    As1   As2   As3   As4  As5  As6   |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*   |   |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*   _   |______________________________________|            *\n');
fprintf(fid,'*       |<-                 h                ->|            *\n');
fprintf(fid,'*       |<-d1->|                                            *\n');
fprintf(fid,'*       |<-  d2   ->|                                       *\n');
fprintf(fid,'*       |<-     d3      ->|                                 *\n');
fprintf(fid,'*       |<-        d4          ->|                          *\n');
fprintf(fid,'*       |<-            d5          ->|                      *\n');
fprintf(fid,'*       |<-               d6             >|                 *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'+==========================+\n');
fprintf(fid,' =       Unconfined       =\n');
fprintf(fid,'  Curvature    Moment      \n');
fprintf(fid,'    (1/m)      (kN.m)      \n');
fprintf(fid,'----------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur;Mom]);
fprintf(fid,'+==========================+\n');
fprintf(fid,' = Unconfined curve fitted =\n');
fprintf(fid,'  Curvature    Moment       \n');
fprintf(fid,'    (1/m)      (kN.m)       \n');
fprintf(fid,'----------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[X;Y]);
fprintf(fid,'+==========================+\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================================================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I;ES;XX;CUR;EI]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);
