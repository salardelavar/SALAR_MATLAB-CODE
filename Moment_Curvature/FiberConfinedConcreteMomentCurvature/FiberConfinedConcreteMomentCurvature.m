%***********************************************************%
%                >> IN THE NAME OF ALLAH <<                 %
%  Moment-Curvature Analysis of Confined Concrete Section   %
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
%       |     #     #     #     #    #    #    |            %
%       |     #                           #    |            %
%   b   |    As1   As2   As3   As4  As5  As6   |            %
%       |     #                           #    |            %
%   |   |     #     #     #     #    #    #    |            %
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
b=400;% [mm]
h=400;% [mm]
%As:  As1      As2     As3     As4    As5      As6
As=[2454.296 0 0 0 0 2454.296]; % NOTE: As1 & As2 = 5fi25
%d:d1  d2  d3  d4  d5  d6 
d=[70.5 0 0 0 0 329.5]; % Distance of rebars - NOTE: d(1)=cover+stirrup diameter+half longitudinal rebar
cover=50; % Concrete cover - NOTE: Cover must less than d(1)
%% Concrete Properties
fc =30;% [N/mm^2] Unconfined Nominal concrete strength
ecu=0.004;% Ultimate unconfined concrete strain
esp=ecu+.002; % Spalling concrete strain
Ec=5000*sqrt(fc);
ec0=(2*fc)/Ec;
fct=-0.7*sqrt(fc);% Concrete crack stress
ect1=(2*fct)/Ec;ect2=(2.625*fct)/Ec;ect3=(9.292*fct)/Ec;%Concrete crack strain
%% Reinforcing steel Properties
fy =400;% [N/mm^2] Yield strength of reinforcing steel
Es =2e5;% [N/mm^2] Modulus of elasticity of steel
fu=1.5*fy;% Ultimate steel stress
ey=fy/Es;% Yield steel strain
esh=0.01;% Strain at steel strain-hardening
esu=0.09;% Ultimate steel strain
stp= 2;% Kind of strain hardening:  1:linear   2:curve
Esh=(fu-fy)/(esu-esh);
%% Reinforcing stirrup Properties
fyh =300;% [MPa] Transverse Reinforcing Bar (stirrup) Yield Stress
Nx=2;% Total of transverse hoop legs in the X
Ny=2;% Total of transverse hoop legs in the Y
diastirrup=8;% [mm] Cross-sectional hoop diameter
wi=360000;% [mm^2] Define total summation of confined reinforcing steel in power.2
s =150;% [mm] Tie Spacing Along Member

N=1000;% Number of concrete Fiber
itermax = 3000;% maximum number of iterations
tolerance = 1e-3; % specified tolerance for convergence
x=.5*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% Calculating concrete properties
Asp=3.1415*(diastirrup^2)/4;
cx=h-2*cover-diastirrup;%confined core dimensions in the X 
cy=b-2*cover-diastirrup;%confined core dimensions in the Y
ro_x=Nx*Asp/(s*cx);ro_y=Ny*Asp/(s*cy);%X and Y transverse steel reinforcing ratio
ro_s=ro_x+ro_y;sBAR=s-diastirrup;
Wi=1-(wi/(6*cx*cy));Ki=(1-(sBAR/(2*cx)))*(1-(sBAR/(2*cy)));
rcc=(sum(As))/(cx*cy);% rcc is the ratio of the longitudinal reinforcement area to section core area
%ke=(Wi*Ki)/(1-rcc);%confinement effective coefficient
ke=0.75;
f_BAR_lx=ke*ro_x*fyh;f_BAR_ly=ke*ro_y*fyh;
f_BAR_lx_fc=f_BAR_lx/fc;
f_BAR_ly_fc=f_BAR_ly/fc;
flx=f_BAR_lx/fc;fly=f_BAR_ly/fc;
% For rectangular concrete section con?ned by rectangular hoops. 
% Peak stress of confined concrete: Mander Solution.
A=196.5*fly^2+29.1*fly-4;
B=-69.5*fly^2+8.9*fly+2.2;
C=-6.83*fly^2+6.38*fly+1;
D=-1.5*fly^2-0.55*fly+0.3;
if (fly < flx) && (fly <= 0.15)
    FCC=A*flx^2+B*fly+C;
elseif (fly < flx) && (fly > 0.15)
    FCC=((flx-fly)*B/(0.3-fly))+C;
else
    FCC=C;
end
fcc=fc*FCC;K=fcc/fc;
ecc=ec0*(1+5*(K-1));
Esec=fcc/ecc;r=Ec/(Ec-Esec);
ecuc=0.004+(1.4*ro_s*fyh*esu)/fcc;% Ultimate confined concrete strain [Empirical Conservative Formulation (Priestley)]
EsecU=fc/ec0;R=Ec/(Ec-EsecU);
fcu=(fc*R*(ecu/ec0))/(R-1+(ecu/ec0)^R);
%% ------------------ Newton Method Procedure ------------------------%
An=size(As,2);T=(1/N);
for k=1:N;c(k)=(.5*T+T*(k-1))*h;end
% EC=[abs(ect1) abs(ect2) .33*abs(ect3) .67*abs(ect3) .8*abs(ect3) .9*abs(ect3) abs(ect3) .4*ecuc .5*ecuc .6*ecuc .7*ecuc .8*ecuc .9*ecuc ecuc];q=size(EC,2);
EC=0.00001:0.0001:ecuc;q=size(EC,2);
for j=1:q;
    ECC=EC(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
    eC=(ECC*x)/(x-cover);% top Unconfined compression strain
             for u=1:An;
             es=eC*(x-d(u))/x;
    %---------------- As -------------%
       if and(es>0,es<ey)
        fs=Es*es;
        fstan=(Es*eC*d(u))/(x)^2;
       elseif and(es<0,es>(-ey))
        fs=Es*es;
        fstan=(Es*eC*d(u))/(x)^2;
       elseif  and(es>=ey,es<esh)
        fs=fy;
        fstan=0;
       elseif and(es<=(-ey),es>(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>=esh,es<esu)
            if stp==1;
               fs=fy+Esh*(abs(es)-esh);
               fstan=(Esh*eC*d(u))/(x)^2;
            elseif stp==2;
               fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
               fstan=(2*eC*d(u)*(fu-fy)*(((eC*d(u))/x)+esu-eC))/(x^2*(esu-esh)^2);
            end
       elseif and(es<=(-esh),es>(-esu))
           if stp==1;
              fs=-fy-Esh*(abs(es)-esh);
              fstan=(Esh*eC*d(u))/(x)^2;
           elseif stp==2;
              fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
              fstan=(2*eC*d(u)*(fu-fy)*(((eC*d(u))/x)+esu-eC))/(x^2*(esu-esh)^2);
           end
       elseif or(es>=esu,es<=(-esu))
        fs=0;
        fstan=0;
     end
        Fs(u)=As(u)*fs;
        Fstan(u)=As(u)*fstan;% tangent steel force
        if  and(es>0,es<ecuc)% in this step: confined concrete force in rebar area is omitted (F=As*fc)
        C_cons=(fcc*r*(es/ecc))/(r-1+(es/ecc)^r);
         c1s=-(d(u)*fcc*r^2*(x-d(u))*es^2*(es*((x-d(u))/(ecc*x)))^(r-1))/(ecc^2*x^3*(((es*(x-d(u)))/(ecc*x))^r +r-1)^2);
         c2s=-(es*fcc*r*(x-d(u)))/(ecc*x^2*(((es*(x-d(u)))/(ecc*x))^r +r-1));
         c3s=(es*fcc*r)/(ecc*x*(((es*(x-d(u)))/(ecc*x))^r +r-1));
         Ctan_cons=c1s+c2s+c3s;
      elseif es>ecuc
         C_cons=0;
         Ctan_cons=0;
       elseif and(es<0,es>=ect1)
        C_cons=-(0.5*Ec*es);
        Ctan_cons=-((0.5*Ec*d(u)*eC)/x^2);
      elseif and(es<ect1,es>=ect2)
        C_cons=-(fct-(0.5*fct/(ect2-ect1))*(es-ect1));
        Ctan_cons=+(0.5*fct*eC*d(u))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        C_cons=-(.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2));
        Ctan_cons=+(0.5*fct*eC*d(u))/((ect3-ect2)*x^2);
      elseif es<ect3
        C_cons=0;
        Ctan_cons=0;  
        end
      C_conS(u)=C_cons;Ctan_conS(u)=Ctan_cons;
        Fss(u)=-As(u)*C_conS(u);
        Fstans(u)=-As(u)*Ctan_conS(u);% tangent Minus of concrete force
             end
    for z=1:N;% in this step: concrete force for each fiber is calculated
     %-------------- Cc --------------%
     ec=eC*(x-c(z))/x;
     if and(z>=1,z<=N*(cover/h)) % TOP SECTIONS
     if and(ec>0,ec<ecu) % unconfined
         C=(fc*R*(ec/ec0))/(R-1+(ec/ec0)^R);
         C1=-(c(z)*fc*R^2*(x-c(z))*ec^2*(ec*((x-c(z))/(ec0*x)))^(R-1))/(ec0^2*x^3*(((ec*(x-c(z)))/(ec0*x))^R +R-1)^2);
         C2=-(ec*fc*R*(x-c(z)))/(ec0*x^2*(((ec*(x-c(z)))/(ec0*x))^R +R-1));
         C3=(ec*fc*R)/(ec0*x*(((ec*(x-c(z)))/(ec0*x))^R +R-1));
         Ctan=C1+C2+C3;
         CD=0;
      elseif and(ec>=ecu,ec<esp)
        C=fcu*(1-((ec-ecu)/(esp-ecu)));
        Ctan=-(eC*c(z)*fcu)/((esp-ecu)*x^2);
        CD=0;
      elseif ec >= esp
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
     Cc(z)=b*T*h*C;% unconfined top 
     Cctan(z)=b*T*h*Ctan;% unconfined top tangent concrete force
     elseif and(z>(cover/h)*N,z<=N*(1-(cover/h))) % MIDDLE SECTIONS
         if and(ec>0,ec<ecu) % unconfined
         C=(fc*R*(ec/ec0))/(R-1+(ec/ec0)^R);
         C1=-(c(z)*fc*R^2*(x-c(z))*ec^2*(ec*((x-c(z))/(ec0*x)))^(R-1))/(ec0^2*x^3*(((ec*(x-c(z)))/(ec0*x))^R +R-1)^2);
         C2=-(ec*fc*R*(x-c(z)))/(ec0*x^2*(((ec*(x-c(z)))/(ec0*x))^R +R-1));
         C3=(ec*fc*R)/(ec0*x*(((ec*(x-c(z)))/(ec0*x))^R +R-1));
         Ctan=C1+C2+C3;
         CD=0;
      elseif and(ec>=ecu,ec<esp)
        C=fcu*(1-((ec-ecu)/(esp-ecu)));
        Ctan=-(eC*c(z)*fcu)/((esp-ecu)*x^2);
        CD=0;
      elseif ec >= esp
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
      if  and(ec>0,ec<ecuc)% confined
        C_con=(fcc*r*(ec/ecc))/(r-1+(ec/ecc)^r);
         c1=-(c(z)*fcc*r^2*(x-c(z))*ec^2*(ec*((x-c(z))/(ecc*x)))^(r-1))/(ecc^2*x^3*(((ec*(x-c(z)))/(ecc*x))^r +r-1)^2);
         c2=-(ec*fcc*r*(x-c(z)))/(ecc*x^2*(((ec*(x-c(z)))/(ecc*x))^r +r-1));
         c3=(ec*fcc*r)/(ecc*x*(((ec*(x-c(z)))/(ecc*x))^r +r-1));
         Ctan_con=c1+c2+c3;
         CD=0;
      elseif ec>ecuc
         C_con=0;
         Ctan_con=0;
         CD=0;
       elseif and(ec<0,ec>=ect1)
        C_con=0.5*Ec*ec;
        Ctan_con=(0.5*Ec*c(z)*eC)/x^2;
        CD=0;
      elseif and(ec<ect1,ec>=ect2)
        C_con=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan_con=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
        CD=0;
      elseif and(ec<ect2,ec>=ect3)
        C_con=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan_con=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
        CD=0;
      elseif ec<ect3
        C_con=0;
        Ctan_con=0;
        CD(z)=h-c(z);% Crack Depth
      end
     Cc(z)=(b-(b-2*cover))*T*h*C+(b-2*cover)*T*h*C_con; %confined and unconfined middle concrete force
     Cctan(z)=(b-(b-2*cover))*T*h*Ctan+(b-2*cover)*T*h*Ctan_con;%confined and unconfined middle tangent concrete force
     elseif and(z>(1-(cover/h))*N,z<=N) % BELOW SECTIONS
         if and(ec>0,ec<ecu) % unconfined
         C=(fc*R*(ec/ec0))/(R-1+(ec/ec0)^R);
         C1=-(c(z)*fc*R^2*(x-c(z))*ec^2*(ec*((x-c(z))/(ec0*x)))^(R-1))/(ec0^2*x^3*(((ec*(x-c(z)))/(ec0*x))^R +R-1)^2);
         C2=-(ec*fc*R*(x-c(z)))/(ec0*x^2*(((ec*(x-c(z)))/(ec0*x))^R +R-1));
         C3=(ec*fc*R)/(ec0*x*(((ec*(x-c(z)))/(ec0*x))^R +R-1));
         Ctan=C1+C2+C3;
         CD=0;
      elseif and(ec>=ecu,ec<esp)
        C=fcu*(1-((ec-ecu)/(esp-ecu)));
        Ctan=-(eC*c(z)*fcu)/((esp-ecu)*x^2);
        CD=0;
      elseif ec >= esp
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
     Cc(z)=b*T*h*C;% unconfined below tangent concrete force
     Cctan(z)=b*T*h*Ctan;% unconfined below tangent concrete force
     end
     Ss(z)=ec;SS(j,z)=Ss(z);% Concrete Fiber Strain
     CFS(j,z)=Cc(z)/(b*T*h);%Concrete Fiber Stress
    end
    %----------------------------------%
    FsTOTAL=sum(Fs);CcTOTAL=sum(Cc);FssTOTAL=sum(Fss);A=CcTOTAL+FsTOTAL+FssTOTAL;
    FsTOTAL_tan=sum(Fstan);CcTOTAL_tan=sum(Cctan);FssTOTAL_tan=sum(Fstans);A_tan=CcTOTAL_tan+FsTOTAL_tan+FssTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eC,A)
          disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    %if it == itermax;break;end % stop the analysis at all because last Convergence     
     for u=1:An;e(u)=x-d(u);end % distance of each rebar from neutral axis
     for k=1:N;cc(k)=x-c(k);end %distance of each concrete fiber from neutral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.6f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eC,x,(eC/x)*1000,(Fs*e'+Cc*cc')*.000001)
        end
BCTS(j)=(ECC*(h-cover-x))/(x-cover);% Bottom Confined tension strain
TCCS(j)=(fcc*r*(ECC/ecc))/(r-1+(ECC/ecc)^r);% Top Confined compression stress
if and(eC>0,eC<ecu);ECCC(j)=eC;TUCS(j)=(fc*R*(eC/ec0))/(R-1+(eC/ec0)^R);elseif and(eC>=ecu,eC<esp);ECCC(j)=eC;TUCS(j)=fcu*(1-((eC-ecu)/(esp-ecu)));elseif eC>=esp;ECCC(j)=eC;TUCS(j)=0;end% Top Unconfined compression stress 
if CD==0;CrackDepth(j)=0; else CrackDepth(j)=max(abs(CD));end;% Crack Depth of each increment
% Calculate Moment and Curavture
Cur(j)=(eC/x)*1000;XX(j)=x;AA(j)=A;I(j)=j;
Mom(j)=roundn((Fs*e'+Fss*e'+Cc*cc'),-3)*.000001;
CUR(j)=Cur(j);EI(j)=Mom(j)/Cur(j);
end
Cur=abs([0 Cur]);Mom=abs([0 Mom]);
s=size(Cur,2);for i=1:s-1;EI(i)=abs((Mom(i))/(Cur(i)));end % Flextural Rigidity
if eC == ecuc;fprintf('\n      ## Confined Concrete Strain Reached to Ultimate Strain: %1.4f ## \n\n',eC);end
%% Confined bilinear fitting
SIZE=size(Mom,2);
for i=1:SIZE-1;
    hh(i) = Cur(i+1)-Cur(i);
    Aa(i)=(Mom(i)+Mom(i+1))*0.5*hh(i);
end
Area=sum(Aa);k0 =Mom(5)/Cur(5);
fiy = (Mom(i+1)*max(Cur)*0.5-Area)/(Mom(i+1)*0.5 - k0*max(Cur)*0.5);
My = k0*fiy;
X = [0 fiy max(Cur)];Y = [0 My Mom(i+1)];
disp('+==========================+');
disp(' =  Confined curve fitted =');
disp('  Curvature    Moment');
disp('    (1/m)      (kN.m)   ');
disp('----------------------------');
disp([X' Y']);
disp('+==========================+');
%% Ductility_Rito of Confined Section
Elastic_EI=Y(2)/X(2);
Plastic_EI=(Y(3)-Y(2))/(X(3)-X(2));
Ductility_Rito=X(3)/X(2);
fprintf('+-----------------------------------------+\n')
fprintf(' Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI)
fprintf(' Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI)
fprintf(' Confined Section Ductility Ratio is: %5.2f\n',Ductility_Rito)
fprintf('+-----------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% XTRACT
CurXT=[0
0.002031
0.004062
0.006094
0.008125
0.01016
0.01219
0.01422
0.01625
0.01828
0.02031
0.02939
0.03848
0.04756
0.05664
0.06572
0.07481
0.08389
0.09297
0.1021
0.1111
0.1202
0.1293
0.1384
0.1475
0.1565
0.1656
0.1747
0.1838
0.1929
0.202
0.2036];

MomXT=[0
65.61
122.6
179.2
233.9
287
288.7
291.3
292.4
293.2
293.4
295.6
294.2
297.4
302.1
302.8
303.1
306.2
310.2
314.7
318.8
324.1
326.2
330.2
334
337.7
340.7
344.6
347.5
350.1
352.4
352.2];
%% C++ OUTPUT:
Curvature=[0
7.51653e-007
1.56108e-006
2.81062e-006
4.11894e-006
5.487e-006
6.83995e-006
8.08921e-006
9.26896e-006
1.19074e-005
1.79062e-005
2.39543e-005
2.95693e-005
3.47554e-005
3.96074e-005
4.35795e-005
4.72596e-005
5.07924e-005
5.40228e-005
5.69965e-005
5.9776e-005
6.23927e-005
6.48724e-005
6.7248e-005
6.95167e-005
7.1711e-005
7.38423e-005
7.59996e-005
7.82145e-005
8.04796e-005
8.28181e-005
8.5226e-005
8.76821e-005
9.01858e-005
9.27355e-005
9.5329e-005
9.79579e-005
0.000100635
0.000103353
0.00010611
0.000108901
0.000111727
0.000114586
0.000117475
0.000120395
0.000123379
0.000126452
0.000129562
0.000132709
0.000135892
0.000139107
0.000142356
0.000145635
0.000148934
0.00015227
0.000155629
0.000159018
0.000162444
0.000165963
0.000169512
0.000173082
0.000176676
0.00017922
0.000180501
0.000181855
0.000183222
0.000184623
0.000185982
0.000187341
0.000188702
0.000190032
0.000191379
0.000192669
0.000194002
0.000195276
0.000196571
0.000197831
0.000199103
0.000200339
0.000201595
0.000202805
0.000204046
0.000205232
0.000206485
0.000207767
0.000209075
0.000210335
0.000211625
0.000212862
0.000214132
0.000215353
0.000216602
0.00021781
0.00021904
0.000220269
0.00022144
0.000222652
0.000223807
0.000225004
0.000226139
0.00022732
];

Moment=[0
4.13415e+007
7.98255e+007
1.08594e+008
1.41291e+008
1.73123e+008
2.04445e+008
2.35184e+008
2.64592e+008
2.822e+008
2.82297e+008
2.82833e+008
2.83114e+008
2.83117e+008
2.8292e+008
2.85037e+008
2.87489e+008
2.89748e+008
2.91491e+008
2.92798e+008
2.93785e+008
2.94527e+008
2.95072e+008
2.95457e+008
2.95712e+008
2.95855e+008
2.95924e+008
2.96114e+008
2.96458e+008
2.96948e+008
2.97438e+008
2.97991e+008
2.98659e+008
2.99436e+008
3.00313e+008
3.01283e+008
3.0234e+008
3.03479e+008
3.04694e+008
3.0598e+008
3.0733e+008
3.08741e+008
3.10207e+008
3.11725e+008
3.13291e+008
3.148e+008
3.16226e+008
3.17701e+008
3.19219e+008
3.20779e+008
3.22377e+008
3.24012e+008
3.2568e+008
3.27378e+008
3.29107e+008
3.30864e+008
3.32647e+008
3.34401e+008
3.35933e+008
3.37485e+008
3.39055e+008
3.40644e+008
3.41701e+008
3.42115e+008
3.42568e+008
3.43033e+008
3.43517e+008
3.43982e+008
3.44445e+008
3.44909e+008
3.45358e+008
3.45814e+008
3.46243e+008
3.46691e+008
3.47112e+008
3.47541e+008
3.47954e+008
3.48371e+008
3.48772e+008
3.49181e+008
3.49569e+008
3.49969e+008
3.50345e+008
3.50721e+008
3.51041e+008
3.5137e+008
3.51681e+008
3.52002e+008
3.52305e+008
3.52618e+008
3.52913e+008
3.53218e+008
3.53508e+008
3.53806e+008
3.54102e+008
3.54378e+008
3.54668e+008
3.54937e+008
3.55221e+008
3.55482e+008
3.5576e+008
];
%% imaging
figure(1)
IMAGE1=imread('FiberConfinedConcreteMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE2=imread('FiberConfinedConcreteMomentCurvature-image2.jpg');
image(IMAGE2);axis image;axis off;
%% Plot
figure(3)
plot(EC,XX,'--','LineWidth',2);xlabel('Top confined concrete fiber strain');ylabel('Neutral axis (mm)')
title(' Plotting of the Top confined concrete strain and neutral axis ','Color','b');grid on
figure(4)
plot(EC,TCCS,ECCC,TUCS,'g--','LineWidth',3);xlabel('Top fiber concrete strain ');ylabel('Top fiber concrete stress (N/mm^2)')
title(' Plotting of the Top fiber confined and unconfined concrete strain-stress during the analysis','Color','b');
legend('Confined','Unconfined','Location','NorthEast');grid on
figure(5)
plot(CrackDepth,EC,'--','LineWidth',2);xlabel('Crack depth (mm)');ylabel('Top concrete strain')
title(' Plotting of the Top concrete strain and Crack depth ','Color','b');grid on
figure(6)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(7)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(8)
plot(CUR,EI,'black--','LineWidth',3)
title('# Flextural Rigidity(EI)-Curvature diagram #','Color','b');
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(9)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral axis-Flextural Rigidity(EI) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(10)
plot(Cur,Mom,X,Y,'r--','LineWidth',3)
title(['# CONFINED SECTION MOMENT CURVATURE ANALYSIS DIAGRAM #',' -  EI : ',int2str(Y(2)/X(2)),' (kN.m^2)  -  EI_t : ',int2str((Y(3)-Y(2))/(X(3)-X(2))),' (kN.m^2)  -  Ductility Rito : ',int2str(Ductility_Rito)],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Confined','Confined bilinear fitted','Location','NorthEastOutside');grid on
figure(11)
plot(Cur*.001,Mom*10^6,'b',CurXT*.001,MomXT*10^6,'r--',Curvature,Moment,'g-.','LineWidth',3)
title(['# CONFINED SECTION MOMENT CURVATURE ANALYSIS DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE');ylabel('MOMENT')
legend('Confined','XTRACT','C++','Location','NorthEastOutside');grid on
%% Output Data to .txt file
fid = fopen('FiberConfinedConcreteMomentCurvature-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*  Moment-Curvature Analysis of Confined Concrete Section   *\n');
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
fprintf(fid,' =        Confined        =\n');
fprintf(fid,'  Curvature    Moment\n');
fprintf(fid,'    (1/m)      (kN.m)   \n');
fprintf(fid,'----------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur;Mom]);
fprintf(fid,'+==========================+\n');
fprintf(fid,' = Confined curve fitted  =\n');
fprintf(fid,'  Curvature    Moment\n');
fprintf(fid,'    (1/m)      (kN.m)   \n');
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
fprintf(fid,'================================================================================\n');
fclose(fid);
