%***********************************************************%
%                >> IN THE NAME OF ALLAH <<                 %
%  Moment-Curvature analysis of FRP with                    %
% Unconfined concrete section                               %
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
%   _    ______________________________________             %
%   |   |                                      |            %
%       |     #     #     #     #    #    #    |||          %
%       |     #                           #    |||          %
%   b   |    As1   As2   As3   As4  As5  As6   |||(FRP)     %
%       |     #                           #    |||          %
%   |   |     #     #     #     #    #    #    |||          %
%   _   |______________________________________|            %
%       |<-                 h                ->|            %
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
%As:  As1      As2     As3     As4    As5      As6
As=[2454.296 0 0 0 0 2454.296]; % NOTE: As1 & As6 = 5fi25
%d:d1  d2  d3  d4  d5  d6 
d=[50 0 0 0 0 450];
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
%% Carbon fber [CFRP] Properties
Ecfrp=62000;% [MPa] CFRP Modulus of elasticity
Fcfrp=930;% [MPa] Yield strength of CFRP
ecfrp=0.015;% Ultimate CFRP strain
Bcfrp=300; %[mm]
Tcfrp=1; %[mm]
%% Glass fber [GFRP] Properties
Egfrp=21000;% [MPa] GFRP Modulus of elasticity
Fgfrp=630;% [MPa] Yield strength of GFRP
egfrp=0.03;% Ultimate GFRP strain
Bgfrp=500; %[mm]
Tgfrp=1.3; %[mm]

N=1000;% Number of concrete Fiber
itermax = 4000;% maximum number of iterations
tolerance = 10e-12; % specified tolerance for convergence
x=.5*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure ------------------------%
N1=9*N/10;% Number of steel section top flange Fiber  
N2=1*N/10;% Number of steel section top flange Fiber   
An=size(As,2);R1=(1/N1);R2=(1/N2);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*h;end % distance of each top concrete section Fiber
for k=1:N2;c2(k)=h+(.5*R2+R2*(k-1))*Tcfrp;end % distance of each FRP section Fiber 
An=size(As,2);c=[c1 c2];
EC=[.15*abs(ect1) .33*abs(ect1) .67*abs(ect1) abs(ect1) abs(ect2) .33*abs(ect3) .67*abs(ect3) .8*abs(ect3) .9*abs(ect3) abs(ect3) .4*ecu .5*ecu .6*ecu .7*ecu .75*ecu .8*ecu .85*ecu .9*ecu .95*ecu ecu];q=size(EC,2);
%EC=10^-4:(10^-4/ecu)*10^-3:ecu;q=size(EC,2);
for j=1:q;
    eC=EC(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
       for u=1:An;
       es=eC*(x-d(u))/x;
    %---------------- As -------------%
       if and(es>=0,es<=ey)
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
        if and(es>=0,es<ec0)% in this step: Unconfined concrete force in rebar area is omitted (F=As*fc)
        Cs=fc*((2*es/ec0)-(es/ec0)^2);
        Ctans=((2*fc)/(ec0^2*x^3))*(ec0*eC*d(u)*x-d(u)*eC^2*x+2*eC^2*d(u)^2);
      elseif and(es>=ec0,es<ecu)
        Cs=fc*(1-(0.15*(es-ec0)/(ecu-ec0)));
        Ctans=-(3*eC*d(u)*fc)/(20*(ecu-ec0)*x^2);
      elseif es >= ecu
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
    for z=1:N;% in this step: Fiber force for each fiber is calculated
    if and(z>=1,z<=N1);
     %-------------- Cc --------------%
     ec=eC*(x-c(z))/x;
     if and(ec>=0,ec<ec0)
        C=fc*((2*ec/ec0)-(ec/ec0)^2);
        Ctan=((2*fc)/(ec0^2*x^3))*(ec0*eC*c(z)*x-c(z)*eC^2*x+2*eC^2*c(z)^2);
      elseif and(ec>=ec0,ec<=ecu)
        C=fc*(1-(0.15*(ec-ec0)/(ecu-ec0)));
        Ctan=-(3*eC*c(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif ec > ecu
         C=0; 
         Ctan=0;
      elseif and(ec<0,ec>=ect1)
        C=0.5*Ec*ec;
        Ctan=(0.5*Ec*c(z)*eC)/x^2;
      elseif and(ec<ect1,ec>=ect2)
        C=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
      elseif and(ec<ect2,ec>=ect3)
        C=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
      elseif ec<ect3
        C=0;
        Ctan=0;
     end
     Cc(z)=b*R1*h*C;
     Cctan(z)=b*R1*h*Ctan;% tangent concrete force
    elseif and(z>N1,z<=N)
     if and(ec>=0,ec<=ecfrp)
        C=Ecfrp*ec;
        Ctan=(Ecfrp*c(z)*eC)/x^2;
     elseif (ec > ecfrp)
        C=0.0;
        Ctan=0.0;
     elseif and(es<0,es>=(-ecfrp))
        C=Ecfrp*ec;
        Ctan=(Ecfrp*c(z)*eC)/x^2;
     else (ec < (-ecfrp))
        C=0.0;
        Ctan=0.0;
     end
     Cc(z)=Bcfrp*R2*Tcfrp*C;
     Cctan(z)=Bcfrp*R2*Tcfrp*Ctan;% tangent CFRP force
    end
     Ss(z)=ec;SS(j,z)=Ss(z);% section Fiber Strain
     CFS(j,z)=C;% Section Fiber Stress
    end
    %----------------------------------%
    FsTOTAL=sum(Fs);CcTOTAL=sum(Cc);FssTOTAL=sum(Fss);A=CcTOTAL+FsTOTAL;
    FsTOTAL_tan=sum(Fstan);CcTOTAL_tan=sum(Cctan);FssTOTAL_tan=sum(Fstans);A_tan=CcTOTAL_tan+FsTOTAL_tan;
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
% Calculate Moment and Curavture
Cur(j)=(eC/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*e'+Fss*e'+Cc*cc')*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
Cur=[0 Cur];Mom=[0 Mom];
s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));end % Flextural Rigidity
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
fprintf('+-----------------------------------------+\n')
fprintf(' Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI)
fprintf(' Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI)
fprintf(' Unconfined Section Ductility Ratio : %5.2f\n',Ductility_Rito)
fprintf('+-----------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% unconfined without FRP
Cur01=[         0
    0.0002
    0.0004
    0.0009
    0.0014
    0.0020
    0.0024
    0.0058
    0.0074
    0.0096
    0.0119
    0.0175
    0.0255
    0.0331
    0.0404
    0.0477
    0.0549
    0.0621];
CUR01=[    0.0002
    0.0004
    0.0009
    0.0014
    0.0020
    0.0024
    0.0058
    0.0074
    0.0096
    0.0119
    0.0175
    0.0255
    0.0331
    0.0404
    0.0477
    0.0549
    0.0621];
Mom01=[         0
   25.3301
   55.4416
  111.4667
  148.5395
  178.5899
  198.8977
  360.5328
  406.0552
  404.7245
  404.5879
  405.4086
  406.5064
  412.9810
  420.1888
  427.1278
  433.9060
  440.5788];
EI01=1.0e+005*[    1.2981
    1.2909
    1.2793
    0.7524
    0.4907
    0.4424
    0.4799
    0.2846
   -0.0062
   -0.0006
    0.0015
    0.0014
    0.0085
    0.0098
    0.0095
    0.0094
    0.0093];
XX01=[  215.2401
  215.6896
  216.5482
  206.0238
  186.4099
  176.6288
  150.3187
  140.6729
  122.5694
  109.7541
   91.3651
   78.5806
   72.6126
   69.2366
   67.0426
   65.5162
   64.3998];
X01=[         0
    0.0029
    0.0621];
Y01=[         0
  375.0638
  440.5788];
%% unconfined with  CFRP
Cur02=[         0
    0.0002
    0.0004
    0.0009
    0.0013
    0.0019
    0.0024
    0.0057
    0.0070
    0.0089
    0.0109
    0.0157
    0.0222
    0.0286
    0.0345
    0.0375
    0.0404
    0.0513
    0.0549
    0.0585
    0.0621];
CUR02=[    0.0002
    0.0004
    0.0009
    0.0013
    0.0019
    0.0024
    0.0057
    0.0070
    0.0089
    0.0109
    0.0157
    0.0222
    0.0286
    0.0345
    0.0375
    0.0404
    0.0513
    0.0549
    0.0585
    0.0621];
Mom02=[         0
   25.5021
   55.8177
  112.2215
  150.6612
  181.5760
  202.5831
  367.8652
  426.3059
  431.5564
  438.4995
  456.6615
  481.6305
  506.2696
  533.3408
  546.5304
  559.5464
  430.5327
  433.9060
  437.2528
  440.5789];
EI02=1.0e+005*[    1.3128
    1.3055
    1.2938
    0.7934
    0.5146
    0.4681
    0.5015
    0.4438
    0.0275
    0.0347
    0.0379
    0.0383
    0.0388
    0.0453
    0.0450
    0.0448
   -0.1177
    0.0094
    0.0093
    0.0093];
XX02=[  216.2102
  216.6616
  217.5241
  207.8764
  188.6837
  179.1313
  153.1266
  148.4854
  131.2765
  119.1577
  101.8630
   89.9584
   83.9900
   81.0537
   80.0473
   79.2415
   66.2161
   65.5163
   64.9183
   64.3994];
X02=[         0
    0.0038
    0.0621];
Y02=[         0
  494.0798
  440.5789];
%% unconfined with  GFRP
Cur03=[         0
    0.0002
    0.0004
    0.0009
    0.0014
    0.0020
    0.0024
    0.0058
    0.0073
    0.0093
    0.0115
    0.0168
    0.0242
    0.0313
    0.0381
    0.0448
    0.0514
    0.0580];
CUR03=[    0.0002
    0.0004
    0.0009
    0.0014
    0.0020
    0.0024
    0.0058
    0.0073
    0.0093
    0.0115
    0.0168
    0.0242
    0.0313
    0.0381
    0.0448
    0.0514
    0.0580];
Mom03=[         0
   25.3887
   55.5696
  111.7236
  149.2656
  179.6140
  200.1630
  363.0486
  413.2710
  414.3873
  416.9169
  424.3759
  434.7867
  448.2434
  463.2949
  477.8928
  492.1990
  506.3033];
EI03=1.0e+005*[    1.3031
    1.2958
    1.2842
    0.7663
    0.4989
    0.4512
    0.4873
    0.3361
    0.0054
    0.0116
    0.0140
    0.0141
    0.0189
    0.0221
    0.0218
    0.0216
    0.0215];
XX03=[  215.5695
  216.0197
  216.8796
  206.6568
  187.1884
  177.4863
  151.2782
  143.4317
  125.6623
  113.1100
   95.1346
   82.6800
   76.6924
   73.4708
   71.4093
   69.9975
   68.9818];
X03=[         0
    0.0029
    0.0580];
Y03=[         0
  376.5345
  506.3033];
%% SAP2000- Moment-Curvature of C50x50-10fi25-CFRP30x1
CurC=[0
    1.299*10^-3
    6.494*10^-3
    7.273*10^-3
    3.506*10^-2
    3.714*10^-2
    6.104*10^-2];
MomC=[0
    149.31
    410.03
    429.5
    536.69
    417.74
    440.04];
%% SAP2000- Moment-Curvature of C50x50-10fi25-GFRP30x1
CurG=[0
    1.169*10^-3
    6.818*10^-3
    1.5*10^-2
    3.097*10^-2
    6.688*10^-2];
MomG=[0
    139.71
    415.26
    423.5
    448.45
    503.82];
%% imaging
figure(1)
IMAGE1=imread('FiberUnconfinedConcreteMomentCurvatureFRP-1.jpg');image(IMAGE1);axis image;axis off;
figure(2)
IMAGE1=imread('FiberUnconfinedConcreteMomentCurvatureFRP-2.jpg');image(IMAGE1);axis image;axis off;
figure(3)
IMAGE1=imread('FiberUnconfinedConcreteMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(4)
IMAGE2=imread('FiberUnconfinedConcreteMomentCurvature-image2.jpg');
image(IMAGE2);axis image;axis off;
%% Plot
figure(5)
p1=plot(I,DX,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
figure(6)
p1=plot(I,IT,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
figure(7)
plot(EC,XX,'--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top concrete strain and Neuteral axis ','Color','b');grid on
figure(8)
plot(EC,CUR,'--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top concrete strain and Curvature ','Color','b');grid on
figure(9)
plot(ECCC,TUCS,'black','LineWidth',3);xlabel('Top fiber concrete strain ');ylabel('Top fiber concrete stress (N/mm^2)')
title(' Plotting of the Top fiber unconfined concrete strain-stress during the analysis','Color','b');grid on
figure(10)
plot(SS(1,:),((h+Tcfrp)-c),SS(roundn(.25*q,0),:),((h+Tcfrp)-c),SS(roundn(.5*q,0),:),((h+Tcfrp)-c),SS(roundn(.75*q,0),:),((h+Tcfrp)-c),SS(roundn(q,0),:),((h+Tcfrp)-c),'g--','LineWidth',3);xlabel('Strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the fiber strain and Heigh of section  with CFRP - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
ylim([0 (h+Tcfrp)])
figure(11)
plot(CFS(1,:),((h+Tcfrp)-c),CFS(roundn(.25*q,0),:),((h+Tcfrp)-c),CFS(roundn(.5*q,0),:),((h+Tcfrp)-c),CFS(roundn(.75*q,0),:),((h+Tcfrp)-c),CFS(roundn(q,0),:),((h+Tcfrp)-c),'g--','LineWidth',3);xlabel('Stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the fiber stress and Heigh of section with CFRP - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
ylim([0 (h+Tcfrp)])
figure(12)
plot(TSCSs,TSCS,BSCSs,BSCS,'g--','LineWidth',3);xlabel('Steel strain (\epsilon_s)');ylabel('Steel stress (N/mm^2)')
title(' Plotting of the Top  and Bottom steel strain-stress during the analysis','Color','b');
legend('Top steel rebar','Bottom steel rebar','Location','NorthEastOutside');grid on;
figure(13)
plot(CUR,EI,'black--','LineWidth',3)
title('# UNCONFINED CONCRETE AND FRP EI-CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(14)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral Axis-EI diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(15)
plot(Cur,Mom,X,Y,'r--','LineWidth',3)
title(['# UNCONFINED CONCRETE AND FRP MOMENT CURVATURE DIAGRAM #',' -  EI_e_l_a : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Unconfined','Unconfined bilinear fitted','Location','NorthEastOutside');grid on;
figure(16)
plot(CUR01,EI01,CUR02,EI02,'--',CUR03,EI03,'--','LineWidth',3)
title('# UNCONFINED CONCRETE AND FRP EI-CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
legend('unconfined sec.','CFRP unconfined sec.','GFRP unconfined sec.','Location','NorthEastOutside');
figure(17)
plot(XX01,EI01,XX02,EI02,'--',XX03,EI03,'--','LineWidth',3)
title('# Neuteral Axis-EI diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
legend('unconfined sec.','CFRP unconfined sec.','GFRP unconfined sec.','Location','NorthEastOutside');
figure(18)
plot(Cur01,Mom01,Cur02,Mom02,'--',Cur03,Mom03,'--','LineWidth',3)
title(['# MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)');grid on;
legend('unconfined sec.','CFRP unconfined sec.','GFRP unconfined sec.','Location','NorthEastOutside');
figure(19)
plot(X01,Y01,X02,Y02,'--',X03,Y03,'--','LineWidth',3)
title(['# MOMENT-CURVATURE CURVE FITTED DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)');grid on;
legend('unconfined sec.','CFRP unconfined sec.','GFRP unconfined sec.','Location','NorthEastOutside');
figure(20)
plot(Cur02,Mom02,CurC,MomC,'r--','LineWidth',3)
title(['# UNCONFINED CONCRETE SECTION AND CFRP MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Analysis - CFRP','SAP2000 - CFRP','Location','NorthEastOutside');grid on;
figure(21)
plot(Cur03,Mom03,CurG,MomG,'r--','LineWidth',3)
title(['# UNCONFINED CONCRETE SECTION AND GFRP MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Analysis - GFRP','SAP2000 - GFRP','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberUnconfinedConcreteMomentCurvatureFRP-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*Moment-Curvature analysis of FRP with Unconfined           *\n');
fprintf(fid,'* concrete section                                          *\n');
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
fprintf(fid,'*       |     #     #     #     #    #    #    |||          *\n');
fprintf(fid,'*       |     #                           #    |||          *\n');
fprintf(fid,'*   b   |    As1   As2   As3   As4  As5  As6   ||| (FRP)    *\n');
fprintf(fid,'*       |     #                           #    |||          *\n');
fprintf(fid,'*   |   |     #     #     #     #    #    #    |||          *\n');
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
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I;EC;XX;CUR;EI]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);
