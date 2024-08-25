%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
%   Moment-Curvature Analysis of Rectangular steel section  %
% With Axial Load effect                                    %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given: Steel section properties ,                          %
%Calculate: Moment-Curavture                                %
% Newton-Raphson Method : Tangent procedure                 %
%***********************************************************%
%   _    ______________________________________             %
%   |   |                                      |            %
%       |                                      |            %
%       |                                      |            %
%   b   |                                      |            %
%       |                                      |            %
%   |   |                                      |            %
%   _   |______________________________________|            %
%       |<-                 h                ->|            %
%    X                                                      %
%    ^                                                      %
%    |             (Moment - Curvature along X axis)        %
%    |                                                      %
%    +----> Y                                               %
%***********************************************************%
clear all;close all;clc;
Ptarget =+2000000;% [N] Target axial load [+ : Compression]
%% Section Properties
b=100;% [mm]
h=100;% [mm]
%% Steel Section Properties
fy =240;% [N/mm^2] Yield strength of steel section
Es =2e5;% [N/mm^2] Modulus of elasticity of steel section
fu=1.5*fy;% Ultimate steel stress
ey=fy/Es;% Yield steel strain
esh=0.025;% Strain at steel strain-hardening
esu=0.35;% Ultimate steel strain
Esh=(fu-fy)/(esu-esh);
N=1000;% Number of steel section Fiber
itermax = 4000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
x=.25*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure ------------------------%
R=(1/N);
for k=1:N;c(k)=(.5*R+R*(k-1))*h;end
if abs(Ptarget)>0
it = 0; % initialize iteration count
residual = 100; % initialize residual
esi = tolerance; % initialize strain guess
    while (residual > tolerance)
    for z=1:N;% in this step: steel force for each fiber is calculated
    %---------------- Fs -------------%
       if and(esi>0,esi<=ey)
        fsi=Es*esi;
        fstani=Es;
       elseif and(esi<0,esi>=(-ey))
        fsi=Es*esi;
        fstani=-Es;
       elseif  and(esi>ey,esi<=esh)
        fsi=fy;
        fstani=0;
       elseif and(esi<(-ey),esi>=(-esh))
        fsi=-fy;
        fstani=0;
        elseif  and(esi>esh,esi<=esu)
        fsi=fu-(fu-fy)*((esu-abs(esi))/(esu-esh))^2;
        fstani=(fu-fy)*(2*(esu-abs(esi)))/(esu-esh)^2;
       elseif and(esi<(-esh),esi>=(-esu))
        fsi=-fu+(fu-fy)*((esu-abs(esi))/(esu-esh))^2;
        fstani=+(fu-fy)*(2*(abs(esi)-esu))/(esu-esh)^2;
       elseif or(esi>esu,esi<(-esu))
        fsi=0;
        fstani=0;
       end
     Fsi(z)=b*R*h*fsi;
     Fstani(z)=b*R*h*fstani;% tangent steel force    
    end
       FsTOTALi=sum(Fsi);Ai=FsTOTALi-Ptarget;
       FsTOTAL_tani=sum(Fstani);A_tani=FsTOTAL_tani;
        dxi = A_tani^-1 *(-Ai);
        residual = abs(dxi); % evaluate residual
        it = it + 1; % increment iteration count
        esi = esi+dxi; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-) trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',it,esi,residual)
          disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
if it < itermax% iteration control
        fprintf('(+)It is converged in %1.0f iterations - Initial axial strain: %1.6f - Initial axial stress: %1.3f (N/mm^2)\n\n',it,esi,fsi)
end
else;esi=0;end
ES=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q=size(ES,2);
% ES=esi+10^-4:((esi+10^-4)/esu)*esi+10^-4:esu;q=size(ES,2);
for j=1:q;
    eS=ES(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=esi+eS*(x-c(z))/x;
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
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
     Fs(z)=b*R*h*fs;
     Fstan(z)=b*R*h*fstan;% tangent steel force
     Ss(z)=es;SS(j,z)=Ss(z);%steel sction Fiber Strain
     CFS(j,z)=fs;%steel sction Fiber Stress            
     end
    %----------------------------------%
    FsTOTAL=sum(Fs);A=FsTOTAL-Ptarget;
    FsTOTAL_tan=sum(Fstan);A_tan=FsTOTAL_tan;
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
    if it == itermax;break;end % stop the analysis at all because last Convergence
     Pc1=x-.5*h;% Distance from applied force from neutral axis
     for k=1:N;cc(k)=x-c(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.4f - x: %1.2f - Phi: %1.4f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc'+Pc1*Ptarget)*.000001)
        end
if and(eS>0,eS<ey);ESS(j)=eS;TUSS(j)=Es*eS;elseif and(eS>=ey,eS<esh);ESS(j)=eS;TUSS(j)=fy;elseif and(eS>=esh,eS<esu);ESS(j)=eS;TUSS(j)=fu-(fu-fy)*((esu-abs(eS))/(esu-esh))^2;end% Top Steel section compression stress        
% Calculate Moment and Curavture
Cur(j)=(eS/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*cc'+Pc1*Ptarget)*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;EI(j)=abs(Mom(j)/Cur(j));
end
Cur=[0 Cur];Mom=[0 Mom];
s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));EA(i)=(12*EI(i))/h^2;end % Flextural and Axial Rigidity
if roundn(eS,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
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
%% Elastic-Perfect Plastic Steel Section Moment-Curvature Analysis (simple) - No Axial Load effect
Ie1=((b*h^3)/12);
phiY1=ey/(.5*h);phiU1=esu/(.5*h);
My1=(fy*Ie1)/(.5*h); % Yeild Moment
Mu1=fy*[(b*h^2)/4];  % Ultimate Moment
EIe1=(My1/phiY1)*1e-9;% Elastic  Flextural Rigidity
EIp1=((Mu1-My1)/(phiU1-phiY1))*1e-9;% Plastic Flextural Rigidity
phi1=[0;phiY1;phiU1]*1000;Mi1=[0;My1;Mu1]*10^-6;
SD1=phiU1/phiY1;a1=Mu1/My1;
disp('+==========================================================+');
disp('= Steel Section curve fitted  = Steel Section Simple Anal. =');
disp('     Curvature    Moment          Curvature    Moment');
disp('       (1/m)      (kN.m)            (1/m)      (kN.m)   ');
disp('------------------------------------------------------------');
disp([X' Y' phi1 Mi1]);
disp('+==========================================================+');
%% EI and Ductility_Rito of Steel Section
Elastic_EI=Y(2)/X(2);
Plastic_EI=(Y(3)-Y(2))/(X(3)-X(2));
Ductility_Rito=X(3)/X(2);
Over_Strength_Factor=Y(3)/Y(2);
Material_Ductility_Rito=esu/ey;
fprintf('+--------------------------------------------------+\n')
fprintf(' Elastic EI (Exact): %5.2f (kN.m^2)\n',Elastic_EI)
fprintf(' Plastic EI (Exact): %5.2f (kN.m^2)\n',Plastic_EI)
fprintf(' Elastic EI (Simple): %5.2f (kN.m^2)\n',EIe1)
fprintf(' Plastic EI (Simple): %5.2f (kN.m^2)\n',EIp1)
fprintf(' Steel Material Ductility Ratio : %5.2f\n',Material_Ductility_Rito)
fprintf(' Steel Section Ductility Ratio (Exact) : %5.2f\n',Ductility_Rito)
fprintf(' Steel Section Over Strength Factor (Exact) : %5.2f\n',Over_Strength_Factor)
fprintf(' Steel Section Ductility Ratio (Simple) : %5.2f\n',SD1)
fprintf(' Steel Section Over Strength Factor (Simple) : %5.2f\n',a1)
fprintf('+--------------------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% imaging
figure(1)
IMAGE1=imread('FiberSteelSectionMomentCurvature-image.jpg');
image(IMAGE1);axis image;axis off;
%% Plot
figure(2)
P1=plot(ESS,TUSS,[0;ey;esu],[0;fy;fy],'r--');set(P1,'LineWidth',3);
xlabel('Strain');ylabel('Stress (kN/mm^2)')
legend('Exact','Elastic-Perfect Plastic (Simple)','Location','NorthEastOutside');grid on;
title('Steel Section Stress-Strain Diagram','FontSize',12,'color','b')
figure(3)
p1=plot(I,DX,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
figure(4)
p1=plot(I,IT,'b--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
figure(5)
plot(ES,XX,'--','LineWidth',3);xlabel('Top steel strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top steel strain and Neuteral axis ','Color','b');grid on
figure(6)
plot(ES,CUR,'--','LineWidth',2);xlabel('Top steel strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top steel strain and Curvature ','Color','b');grid on
figure(7)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Steel strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the steel fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(8)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Steel stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the steel fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(9)
plot(CUR,EI,'black--','LineWidth',3)
title('#  Flextural Rigidity(EI)-CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(10)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral Axis-Flextural Rigidity(EI) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(11)
plot(XX,EA,'black','LineWidth',3)
title('# Neuteral Axis-Axial Rigidity(EA) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EA (kN)');grid on;
figure(12)
plot(Cur,Mom,X,Y,'r--',phi1,Mi1,'g-.','LineWidth',3)
title(['# STEEL SECTION MOMENT CURVATURE DIAGRAM WITH AXIAL LOAD EFFECT #','     Target axial load : ',int2str(Ptarget*.001),' (kN) EI_e_l_a : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Exact','Exact - bilinear fitted','Elastic-Perfect Plastic (Simple) - No Axial Load','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberRecSteelSectionMomentCurvatureAxialLoad-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*   Moment-Curvature Analysis of Rectangular steel section  *\n');
fprintf(fid,'* With Axial Load effect                                    *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*     This program is written by salar delavar ghashghaei   *\n');  
fprintf(fid,'*          E-mail:salar.d.ghashghaei@gmail.com              *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*Unit: Newton-Milimeter                                     *\n');
fprintf(fid,'*Given:Section Properties , Steel Section properties ,      *\n');
fprintf(fid,'*Calculate: Moment-Curavture                                *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*   _    ______________________________________             *\n');
fprintf(fid,'*   |   |                                      |            *\n');
fprintf(fid,'*       |                                      |            *\n');
fprintf(fid,'*       |                                      |            *\n');
fprintf(fid,'*   b   |                                      |            *\n');
fprintf(fid,'*       |                                      |            *\n');
fprintf(fid,'*   |   |                                      |            *\n');
fprintf(fid,'*   _   |______________________________________|            *\n');
fprintf(fid,'*       |<-                 h                ->|            *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================+\n');
fprintf(fid,'=Steel Section Moment-Curvature=\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur;Mom]);
fprintf(fid,'+=============================+\n');
fprintf(fid,'= Steel Section curve fitted  =\n');
fprintf(fid,'     Curvature    Moment       \n');
fprintf(fid,'       (1/m)      (kN.m)       \n');
fprintf(fid,'-------------------------------\n');
fprintf(fid,'%10.4f %10.3f\n',[X;Y]);
fprintf(fid,'+=============================+\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================================================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I;ES;XX;CUR;EI]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);