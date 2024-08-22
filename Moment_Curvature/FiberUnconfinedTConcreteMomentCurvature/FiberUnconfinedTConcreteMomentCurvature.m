%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
% Moment-Curvature Analysis of Unconfined T concrete section%
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given:Section Properties , Concrete properties ,           %
% Reinforcing steel properties                              %
%Calculate: Moment-Curavture                                %
% Note: No limit for accounting plurality steel rebar       %
% Newton-Raphson Method : Tangent procedure                 %
%***********************************************************%
%        _     __                                           %   
%        ^    |  |                                          %
%        |    |  |                                          %
%             |  |                                          %
%             |  |_______________ +                         %
%       bf    |   _______________|  tw                      %
%             |  |                +                         %
%        |    |  |                                          %
%        v    |  |                                          %
%        -     ---                                          %
%             |<-   tf1  + hw   ->|                         % 
%    X                                                      %
%    ^                                                      %
%    |             (Moment - Curvature along X axis)        %
%    |                                                      %
%    +----> Y                                               %
%***********************************************************%
%        ______________________________________             %
%       |                                      |            %
%       |     #     #     #     #    #    #    |            %
%       |     #                           #    |            %
%       |    As1   As2   As3   As4  As5  As6   |            %
%       |     #                           #    |            %
%       |     #     #     #     #    #    #    |            %
%       |______________________________________|            %
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
tf=50;% [mm] I section thickness on Top flange
bf=500;% [mm] I section width on Top flange
tw=100;% [mm] I section thickness of Web
hw=200;% [mm] Height of web
tf2=0;% [mm] I section thickness on Bottom flange
bf2=0;% [mm] I section width on Bottom flange
h=tf+tf2+hw;% [mm] Height of Section
%As:  As1      As2     As3     As4    As5      As6
As=[78.5 0 0 0 0 157];       % 3fi10 
%d:d1  d2  d3  d4  d5  d6 
d=[25 0 0 0 0 225];
%% Concrete Properties
fc =25;% [N/mm^2] Unconfined concrete strength 
ecu=0.004;% Ultimate concrete strain
Ec=5000*sqrt(fc);
ec0=(2*fc)/Ec;
fct=-0.7*sqrt(fc);% Concrete tension stress
ect1=(2*fct)/Ec;ect2=(2.625*fct)/Ec;ect3=(9.292*fct)/Ec;% Concrete tension strain
%% Reinforcing steel Properties
fy =400;% [N/mm^2] Yield strength of reinforcing steel
Es =2e5;% [N/mm^2] Modulus of elasticity of steel
fu=1.5*fy;% Ultimate steel stress
ey=fy/Es;% Yield steel strain
esh=0.01;% Strain at steel strain-hardening
esu=0.09;% Ultimate steel strain
Esh=(fu-fy)/(esu-esh);
N=h;% Number of concrete Fiber
itermax = 4000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
x=.5*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure ------------------------%
An=size(As,2);
N1=N/10;% Number of steel section top flange Fiber  
N2=8*N/10;% Number of steel section  web Fiber
N3=N/10;% Number of steel section bottom flange Fiber 
R1=(1/N1);R2=(1/N2);R3=(1/N3);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*tf;end % distance of each top flange Fiber 
for k=1:N2;c2(k)=tf+(.5*R2+R2*(k-1))*hw;end % distance of each web Fiber
for k=1:N3;c3(k)=tf+hw+(.5*R3+R3*(k-1))*tf2;end % distance of each bottom flange Fiber
c=[c1 c2 c3];
EC=[.25*abs(ect1) .5*abs(ect1) .75*abs(ect1) abs(ect1) abs(ect2) .33*abs(ect3) .67*abs(ect3) .8*abs(ect3) .9*abs(ect3) abs(ect3) .4*ecu .5*ecu .6*ecu .7*ecu .8*ecu .9*ecu ecu];q=size(EC,2);
% EC=10^-5:(10^-5/ecu)*10^-5:ecu;q=size(EC,2);
for j=1:q;
    eC=EC(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
       for u=1:An;
       es=eC*(x-d(u))/x;
    %---------------- As -------------%
       if and(es>0,es<ey)
        fs(u)=Es*es;
        fstan(u)=(Es*eC*d(u))/(x)^2;
       elseif and(es<0,es>(-ey))
        fs(u)=Es*es;
        fstan(u)=(Es*eC*d(u))/(x)^2;
       elseif  and(es>=ey,es<esh)
        fs(u)=fy;
        fstan(u)=0;
       elseif and(es<=(-ey),es>(-esh))
        fs(u)=-fy;
        fstan(u)=0;
        elseif  and(es>=esh,es<esu)
        fs(u)=fy+Esh*(abs(es)-esh);
        fstan(u)=(Esh*eC*d(u))/(x)^2;
       elseif and(es<=(-esh),es>(-esu))
        fs(u)=-fy-Esh*(abs(es)-esh);
        fstan(u)=(Esh*eC*d(u))/(x)^2;
       elseif or(es>=esu,es<=(-esu))
        fs(u)=0;
        fstan(u)=0;
     end
        Fs(u)=As(u)*fs(u);
        Fstan(u)=As(u)*fstan(u);% tangent steel force
        if and(es>0,es<ec0)% in this step: Unconfined concrete force in rebar area is omitted (F=As*fc)
        Cs=fc*((2*es/ec0)-(es/ec0)^2);
        Ctans=((2*fc)/(ec0^2*x^3))*(ec0*eC*d(u)*x-d(u)*eC^2*x+2*eC^2*d(u)^2);
      elseif and(es>=ec0,es<=ecu)
        Cs=fc*(1-(0.15*(es-ec0)/(ecu-ec0)));
        Ctans=-(3*eC*d(u)*fc)/(20*(ecu-ec0)*x^2);
      elseif es > ecu
         Cs=0; 
         Ctans=0;
      elseif and(es<0,es>=ect1)
        Cs=-0.5*Ec*es;
        Ctans=-(0.5*Ec*d(u)*eC)/x^2;
      elseif and(es<ect1,es>=ect2)
        Cs=fct+(0.5*fct/(ect2-ect1))*(es-ect1);
        Ctans=+(0.5*fct*eC*d(u))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        Cs=-(.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2));
        Ctans=+(0.5*fct*eC*d(u))/((ect3-ect2)*x^2);
      elseif es<ect3
        Cs=0;
        Ctans=0;
        end
        CS(u)=Cs;CtanS(u)=Ctans;
        Fss(u)=-As(u)*CS(u);
        Fstans(u)=-As(u)*CtanS(u);% tangent Minus of concrete force
    end                
    for z=1:N;% in this step: concrete force for each fiber is calculated
     %-------------- Cc --------------%
     ec=eC*(x-c(z))/x;
     if and(z>=1,z<=roundn(N/10,0)) % TOP FLANGE
     if and(ec>0,ec<ec0)
        C=fc*((2*ec/ec0)-(ec/ec0)^2);
        Ctan=((2*fc)/(ec0^2*x^3))*(ec0*eC*c(z)*x-c(z)*eC^2*x+2*eC^2*c(z)^2);
        CD=0;
      elseif and(ec>=ec0,ec<=ecu)
        C=fc*(1-(0.15*(ec-ec0)/(ecu-ec0)));
        Ctan=-(3*eC*c(z)*fc)/(20*(ecu-ec0)*x^2);
        CD=0;
      elseif ec > ecu
         C=0; 
         Ctan=0;
         CD=0;
      elseif and(ec<0,ec>=ect1)
        C=0.5*Ec*ec;
        Ctan=(0.5*Ec*c(z)*eC)/x^2;
        CD=0;
      elseif and(ec<ect1,ec>=ect2)
        C=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
        CD=0;
      elseif and(ec<ect2,ec>=ect3)
        C=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
        CD=0;
      elseif ec<ect3
        C=0;
        Ctan=0;
        CD(z)=h-c(z);% Crack Depth
     end
     Cc(z)=bf*R1*tf*C;
     Cctan(z)=bf*R1*tf*Ctan;% tangent concrete force
     elseif and(z>roundn(N/10,0),z<=roundn(9*N/10,0)) % MIDDLE WEB
     if and(ec>0,ec<ec0)
        C=fc*((2*ec/ec0)-(ec/ec0)^2);
        Ctan=((2*fc)/(ec0^2*x^3))*(ec0*eC*c(z)*x-c(z)*eC^2*x+2*eC^2*c(z)^2);
        CD=0;
      elseif and(ec>=ec0,ec<=ecu)
        C=fc*(1-(0.15*(ec-ec0)/(ecu-ec0)));
        Ctan=-(3*eC*c(z)*fc)/(20*(ecu-ec0)*x^2);
        CD=0;
      elseif ec > ecu
         C=0; 
         Ctan=0;
         CD=0;
      elseif and(ec<0,ec>=ect1)
        C=0.5*Ec*ec;
        Ctan=(0.5*Ec*c(z)*eC)/x^2;
        CD=0;
      elseif and(ec<ect1,ec>=ect2)
        C=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
        CD=0;
      elseif and(ec<ect2,ec>=ect3)
        C=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
        CD=0;
      elseif ec<ect3
        C=0;
        Ctan=0;
        CD(z)=h-c(z);% Crack Depth
     end
     Cc(z)=tw*R2*hw*C;
     Cctan(z)=tw*R2*hw*Ctan;% tangent concrete force
     elseif and(z>roundn(9*N/10,0),z<=N) % BELOW FLANGE
              if and(ec>0,ec<ec0)
        C=fc*((2*ec/ec0)-(ec/ec0)^2);
        Ctan=((2*fc)/(ec0^2*x^3))*(ec0*eC*c(z)*x-c(z)*eC^2*x+2*eC^2*c(z)^2);
        CD=0;
      elseif and(ec>=ec0,ec<=ecu)
        C=fc*(1-(0.15*(ec-ec0)/(ecu-ec0)));
        Ctan=-(3*eC*c(z)*fc)/(20*(ecu-ec0)*x^2);
        CD=0;
      elseif ec > ecu
         C=0; 
         Ctan=0;
         CD=0;
      elseif and(ec<0,ec>=ect1)
        C=0.5*Ec*ec;
        Ctan=(0.5*Ec*c(z)*eC)/x^2;
        CD=0;
      elseif and(ec<ect1,ec>=ect2)
        C=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
        CD=0;
      elseif and(ec<ect2,ec>=ect3)
        C=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
        CD=0;
      elseif ec<ect3
        C=0;
        Ctan=0;
        CD(z)=h-c(z);% Crack Depth
     end
     Cc(z)=bf2*R3*tf2*C;
     Cctan(z)=bf2*R3*tf2*Ctan;% tangent concrete force
     end
     Ss(z)=ec;SS(j,z)=Ss(z);% Concrete Fiber Strain
     CFS(j,z)=C;% Concrete Fiber Stress
    end                
    %----------------------------------%
    FsTOTAL=sum(Fs);CcTOTAL=sum(Cc);FssTOTAL=sum(Fss);A=CcTOTAL+FsTOTAL+FssTOTAL;
    FsTOTAL_tan=sum(Fstan);CcTOTAL_tan=sum(Cctan);FssTOTAL_tan=sum(Fstans);A_tan=CcTOTAL_tan+FsTOTAL_tan+FssTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = max(abs(dx)); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eC,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence
    for u=1:An;e(u)=x-d(u);end % distance of each rebar from Neuteral axis
     for k=1:N;cc(k)=x-c(k);end %distance of each concrete fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.6f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eC,x,(eC/x)*1000,(Fs*e'+Fss*e'+Cc*cc')*10^-6)
        end
        TSCSs(j)=eC*(x-d(1))/x;TSCS(j)=fs(1);% Top Steel compression strain-stress
        BSCSs(j)=eC*(x-d(An))/x;BSCS(j)=fs(An);% Bottom Steel compression strain-stress
if and(eC>0,eC<ec0);ECCC(j)=eC;TUCS(j)=fc*((2*eC/ec0)-(eC/ec0)^2);elseif and(eC>=ec0,eC<=ecu);ECCC(j)=eC;TUCS(j)=fc*(1-(0.15*(eC-ec0)/(ecu-ec0)));end% Top Unconfined compression stress 
if CD==0;CrackDepth(j)=0; else CrackDepth(j)=max(abs(CD));end;% Crack Depth of each increment
% Calculate Moment and Curavture
Cur(j)=(eC/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*e'+Fss*e'+Cc*cc')*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
Cur=[0 Cur];Mom=[0 Mom];
s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i))/(Cur(i));end % Flextural Rigidity
if roundn(eC,-5) == ecu;fprintf('\n      ## Unconfined Concrete Strain Reached to Ultimate Strain: %1.4f ## \n\n',eC);end
%% UnConfined bilinear fitting
SIZE=size(Mom,2);
for i=1:SIZE-1;
    hh(i) = Cur(i+1)-Cur(i);
    Aa(i)=(Mom(i)+Mom(i+1))*0.5*hh(i);
end
Area=sum(Aa);k0 =Mom(2)/Cur(2);
fiy = (Mom(i+1)*max(Cur)*0.5-Area)/(Mom(i+1)*0.5 - k0*max(Cur)*0.5);
My = k0*fiy;
X = [0 fiy max(Cur)];Y = [0 My Mom(i+1)];
disp('+==========================+');
disp('= Unconfined curve fitted =');
disp('  Curvature    Moment');
disp('    (1/m)      (kN.m)   ');
disp('----------------------------');
disp([X' Y']);
disp('+==========================+');
%% EI and Ductility_Rito of  Unconfined Section
Elastic_EI=Y(2)/X(2);
Plastic_EI=(Y(3)-Y(2))/(X(3)-X(2));
Ductility_Rito=X(3)/X(2);
fprintf('+----------------------------------------+\n')
fprintf(' Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI)
fprintf(' Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI)
fprintf(' Unconfined Section Ductility Ratio : %5.2f\n',Ductility_Rito)
fprintf('+----------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% Moment-Curvature of T25-3fi10
CurP=[0
    2.597*10^-3
    1.169*10^-2
    4.416*10^-2
    5.584*10^-2
    1.714*10^-1
    3.24*10^-1];
MomP=[0
    8.03
    14.42
    14.19
    14.3
    16.8
    19.55];
%% imaging
figure(1)
IMAGE1=imread('FiberUnconfinedConcreteMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE2=imread('FiberUnconfinedTConcreteMomentCurvature.jpg');
image(IMAGE2);axis image;axis off;
%% Plot
figure(3)
p1=plot(I,DX,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
figure(4)
p1=plot(I,IT,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
figure(5)
plot(CrackDepth,EC,'--','LineWidth',2);xlabel('Crack depth (mm)');ylabel('Top concrete strain')
title(' Plotting of the Top concrete strain and Crack depth ','Color','b');grid on
figure(6)
plot(EC,XX,'--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top concrete strain and Neuteral axis ','Color','b');grid on
figure(7)
plot(EC,CUR,'--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top concrete strain and Curvature ','Color','b');grid on
figure(8)
plot(ECCC,TUCS,'black','LineWidth',3);xlabel('Top fiber concrete strain ');ylabel('Top fiber concrete stress (N/mm^2)')
title(' Plotting of the Top fiber unconfined concrete strain-stress during the analysis','Color','b');grid on
figure(9)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(10)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(11)
plot(TSCSs,TSCS,BSCSs,BSCS,'g--','LineWidth',3);xlabel('Steel strain (\epsilon_s)');ylabel('Steel stress (N/mm^2)')
title(' Plotting of the Top  and Bottom steel strain-stress during the analysis','Color','b');
legend('Top steel','Bottom steel','Location','NorthEastOutside');grid on;
figure(12)
plot(CUR,EI,'black','LineWidth',3)
title('# Unconfined Flextural Rigidity(EI)-Curvature Diagram #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(13)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral Axis-Flextural Rigidity(EI) Diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(14)
plot(Cur,Mom,X,Y,'r--','LineWidth',3)
title(['# UNCONFINED MOMENT CURVATURE DIAGRAM #',' -  EI_e_l_a : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Unconfined','Unconfined bilinear fitted','Location','NorthEastOutside');grid on;
figure(15)
plot(Cur,Mom,CurP,MomP,'r--','LineWidth',3)
title(['# UNCONFINED CONCRETE SECTION MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Analysis','SAP2000','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberUnconfinedTConcreteMomentCurvature-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'* Moment-Curvature Analysis of Unconfined T concrete section*\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*     This program is written by salar delavar ghashghaei   *\n');  
fprintf(fid,'*          E-mail:salar.d.ghashghaei@gmail.com              *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*Unit: Newton-Milimeter                                     *\n');
fprintf(fid,'*Given:Section Properties , Concrete properties ,           *\n');
fprintf(fid,'* Reinforcing steel properties                              *\n');
fprintf(fid,'*Calculate: Moment-Curavture                                *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*        _     __                 __                        *\n');   
fprintf(fid,'*        ^    |  |               |  |                       *\n');
fprintf(fid,'*        |    |  |               |  |                       *\n');
fprintf(fid,'*             |  |               |  |                       *\n');
fprintf(fid,'*             |  |_______________|  |  +                    *\n');
fprintf(fid,'*       bf    |   _______________   |  tw                   *\n');
fprintf(fid,'*             |  |               |  |  +                    *\n');
fprintf(fid,'*        |    |  |               |  |                       *\n');
fprintf(fid,'*        v    |  |               |  |                       *\n');
fprintf(fid,'*        -     ---               ---                        *\n');
fprintf(fid,'*             |<-    2*tf+ hw     ->|                       *\n'); 
fprintf(fid,'*    X                                                      *\n');
fprintf(fid,'*    ^                                                      *\n');
fprintf(fid,'*    |             (Moment - Curvature along X axis)        *\n');
fprintf(fid,'*    |                                                      *\n');
fprintf(fid,'*    +----> Y                                               *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*        ______________________________________             *\n');
fprintf(fid,'*       |                                      |            *\n');
fprintf(fid,'*       |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*       |    As1   As2   As3   As4  As5  As6   |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*       |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*       |______________________________________|            *\n');
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
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I;EC;XX;CUR;EI]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);