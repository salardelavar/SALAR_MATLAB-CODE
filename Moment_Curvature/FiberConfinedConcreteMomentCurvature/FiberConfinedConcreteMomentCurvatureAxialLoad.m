%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
%   Moment-Curvature Analysis of Confined concrete section  %
%   With Axial Load effect                                  %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given:Section Properties , Concrete properties ,           %
% Reinforcing steel properties                              %
%Calculate: Moment-Curavture                                %
% Note: No limit for accounting plurality steel rebar       %
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
Ptarget =0;%1000000;% [N] Target axial load [+ : Compression]
%% Section Properties
b=400;% [mm]
h=400;% [mm]
%As:  As1      As2     As3     As4    As5      As6
As=[2454.296 0 0 0 0 2454.296]; % NOTE: As1 & As2 = 5fi25
%d:d1  d2  d3  d4  d5  d6 
d=[70.5 0 0 0 0 329.5]; % Distance of rebars - NOTE: d(1)=cover+stirrup diameter+half longitudinal rebar
cover=50; % Concrete cover - NOTE: Cover must less than d(1)
%% Concrete Properties
fc =30;% [N/mm^2] Unconfined concrete strength 
ecu=0.004;% Ultimate concrete strain
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
fyh =300;%input('Define Transverse Reinforcing Bar (stirrup) Yield Stress (MPa): ');
Nx=2;% input('Define total of transverse hoop legs in the X: ');
Ny=2;% input('Define total of transverse hoop legs in the Y: ');
diastirrup=8;% input('Define cross-sectional hoop diameter: ');
wi=360000;% input('Define total summation of confined reinforcing steel in power.2 (mm^2): ');
s =150;% input('Define Tie Spacing Along Member (mm): ');

N=800;% Number of concrete Fiber
itermax = 500;% maximum number of iterations
tolerance =1e-12; % specified tolerance for convergence
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
if abs(Ptarget)>0
it = 0; % initialize iteration count
residual = 100; % initialize residual
ECCi=1e-12;
    while (residual > tolerance)
             for u=1:An;
    %---------------- As -------------%
       if and(ECCi>=0,ECCi<=ey)
        fsi(u)=Es*ECCi;
        fstani(u)=Es;
       elseif and(ECCi<=0,ECCi>=(-ey))
        fsi(u)=Es*ECCi;
        fstani(u)=-Es;
       elseif  and(ECCi>ey,ECCi<=esh)
        fsi(u)=fy;
        fstani(u)=0;
       elseif and(ECCi<(-ey),ECCi>=(-esh))
        fsi(u)=-fy;
        fstani(u)=0;
        elseif  and(ECCi>esh,ECCi<=esu)
            if stp==1;
               fsi(u)=fy+Esh*(abs(ECCi)-esh);
               fstani(u)=Esh;
            elseif stp==2;
               fsi(u)=fu-(fu-fy)*((esu-abs(ECCi))/(esu-esh))^2;
               fstani(u)=(fu-fy)*(2*(esu-abs(ECCi)))/(esu-esh)^2;
            end
       elseif and(es<(-esh),es>=(-esu))
           if stp==1;
              fsi(u)=-fy-Esh*(abs(ECCi)-esh);
              fstani(u)=-Esh;
           elseif stp==2;
              fsi(u)=-fu+(fu-fy)*((esu-abs(ECCi))/(esu-esh))^2;
              fstani(u)=+(fu-fy)*(2*(abs(ECCi)-esu))/(esu-esh)^2;
           end
       elseif or(es>=esu,es<=(-esu))
        fsi(u)=0;
        fstani(u)=0;
     end
        Fsi(u)=As(u)*fsi(u);
        Fstani(u)=As(u)*fstani(u);% tangent steel force
        if  and(ECCi>0,ECCi<ecuc)% in this step: confined concrete force in rebar area is omitted (F=As*fc)
        C_consi(u)=(fcc*r*(ECCi/ecc))/(r-1+(ECCi/ecc)^r);
         c1si=-(fcc*r^2*ECCi*(ECCi/ecc)^(r-1))/(ecc^2*((ECCi/ecc)^r +r-1)^2);
         c2si=+(fcc*r)/(ecc*((ECCi/ecc)^r +r-1));
         Ctan_consi(u)=c1si+c2si;
      elseif ECCi>ecuc
         C_consi(u)=0;
         Ctan_consi(u)=0;
       elseif and(ECCi<0,ECCi>=ect1)
        C_consi(u)=-(0.5*Ec*ECCi);
        Ctan_consi(u)=-0.5*Ec;
      elseif and(ECCi<ect1,ECCi>=ect2)
        C_consi(u)=-(fct-(0.5*fct/(ect2-ect1))*(ECCi-ect1));
        Ctan_consi(u)=+(0.5*fct*eC*d(u))/((ect2-ect1)*x^2);
      elseif and(ECCi<ect2,ECCi>=ect3)
        C_consi(u)=-(.5*fct-(0.5*fct/(ect3-ect2))*(ECCi-ect2));
        Ctan_consi(u)=+(0.5*fct*eC*d(u))/((ect3-ect2)*x^2);
      elseif ECCi<ect3
        C_consi(u)=0;
        Ctan_consi(u)=0;  
        end
        Fssi(u)=-As(u)*C_consi(u);
        Fstansi(u)=-As(u)*Ctan_consi(u);% tangent Minus of concrete force
             end
    for z=1:N;% in this step: concrete force for each fiber is calculated
     %-------------- Cc --------------%
     if and(z>=1,z<=roundn(N*(cover/h),0)) % TOP SECTIONS
     if and(ECCi>0,ECCi<ecu) % unconfined
         Ci=(fc*R*(ECCi/ec0))/(R-1+(ECCi/ec0)^R);
         C1i=-(fc*R^2*ECCi*(ECCi/ec0)^(R-1))/(ec0^2*((ECCi/ec0)^R +R-1)^2);
         C2i=+(fc*R)/(ec0*((ECCi/ec0)^R +R-1));
         Ctani=C1i+C2i;
      elseif and(ECCi>=ecu,ECCi<esp)
        Ci=fcu*(1-((ECCi-ecu)/(esp-ecu)));
        Ctani=-fcu/(esp-ecu);
      elseif ECCi >= esp
         Ci=0; 
         Ctani=0;
      elseif and(ECCi<0,ECCi>=ect1)
        Ci=0.5*Ec*ECCi;
        Ctani=-0.5*Ec;
      elseif and(ECCi<ect1,ECCi>=ect2)
        Ci=fct-(0.5*fct/(ect2-ect1))*(ECCi-ect1);
        Ctani=-(0.5*fct/(ect2-ect1));
      elseif and(ECCi<ect2,ECCi>=ect3)
        Ci=.5*fct-(0.5*fct/(ect3-ect2))*(ECCi-ect2);
        Ctani=-(0.5*fct/(ect3-ect2));
      elseif ECCi<ect3
        Ci=0;
        Ctani=0;
     end
     Cci(z)=b*T*h*Ci;% unconfined top 
     Cctani(z)=b*T*h*Ctani;% unconfined top tangent concrete force
     elseif and(z>roundn((cover/h)*N,0),z<=roundn(N*(1-(cover/h)),0)) % MIDDLE SECTIONS
         if and(ECCi>0,ECCi<ecu) % unconfined
         Ci=(fc*R*(ECCi/ec0))/(R-1+(ECCi/ec0)^R);
         C1i=-(fc*R^2*ECCi*(ECCi/ec0)^(R-1))/(ec0^2*((ECCi/ec0)^R +R-1)^2);
         C2i=+(fc*R)/(ec0*((ECCi/ec0)^R +R-1));
         Ctani=C1i+C2i;
      elseif and(ECCi>=ecu,ECCi<esp)
        Ci=fcu*(1-((ec-ecu)/(esp-ecu)));
        Ctani=-fcu/(esp-ecu);
      elseif ECCi >= esp
         Ci=0; 
         Ctani=0;
      elseif and(ECCi<0,ECCi>=ect1)
        Ci=0.5*Ec*ECCi;
        Ctani=-0.5*Ec;
      elseif and(ECCi<ect1,ECCi>=ect2)
        Ci=fct-(0.5*fct/(ect2-ect1))*(ECCi-ect1);
        Ctani=-(0.5*fct/(ect2-ect1));
      elseif and(ECCi<ect2,ECCi>=ect3)
        Ci=.5*fct-(0.5*fct/(ect3-ect2))*(ECCi-ect2);
        Ctani=-(0.5*fct/(ect3-ect2));
      elseif ECCi<ect3
        Ci=0;
        Ctani=0;
         end
      if  and(ECCi>0,ECCi<ecuc)% confined
        C_coni=(fcc*r*(ECCi/ecc))/(r-1+(ECCi/ecc)^r);
         c1i=-(fcc*r^2*ECCi*(ECCi/ecc)^(r-1))/(ecc^2*((ECCi/ecc)^r +r-1)^2);
         c2i=+(fcc*r)/(ecc*((ECCi/ecc)^r +r-1));
         Ctan_coni=c1i+c2i;
      elseif ECCi>ecuc
         C_coni=0;
         Ctan_coni=0;
       elseif and(ECCi<0,ECCi>=ect1)
        C_coni=0.5*Ec*ECCi;
        Ctan_coni=-0.5*Ec;
      elseif and(ECCi<ect1,ECCi>=ect2)
        C_coni=fct-(0.5*fct/(ect2-ect1))*(ECCi-ect1);
        Ctan_coni=-(0.5*fct/(ect2-ect1));
      elseif and(ECCi<ect2,ECCi>=ect3)
        C_coni=.5*fct-(0.5*fct/(ect3-ect2))*(ECCi-ect2);
        Ctan_coni=-(0.5*fct/(ect3-ect2));
      elseif ECCi<ect3
        C_coni=0;
        Ctan_coni=0;  
      end
     Cci(z)=(b-(b-2*cover))*T*h*Ci+(b-2*cover)*T*h*C_coni; %confined and unconfined middle concrete force
     Cctani(z)=(b-(b-2*cover))*T*h*Ctani+(b-2*cover)*T*h*Ctan_coni;%confined and unconfined middle tangent concrete force
     elseif and(z>roundn((1-(cover/h))*N,0),z<=N) % BELOW SECTIONS
         if and(ECCi>0,ECCi<ecu) % unconfined
         Ci=(fc*R*(ECCi/ec0))/(R-1+(ECCi/ec0)^R);
         C1i=-(fc*R^2*ECCi*(ECCi/ec0)^(R-1))/(ec0^2*((ECCi/ec0)^R +R-1)^2);
         C2i=+(fc*R)/(ec0*((ECCi/ec0)^R +R-1));
         Ctani=C1i+C2i;
      elseif and(ECCi>=ecu,ECCi<esp)
        Ci=fcu*(1-((ECCi-ecu)/(esp-ecu)));
        Ctani=-fcu/(esp-ecu);
      elseif ECCi >= esp
         Ci=0; 
         Ctani=0;
      elseif and(ECCi<0,ECCi>=ect1)
        Ci=0.5*Ec*ECCi;
        Ctani=-0.5*Ec;
      elseif and(ECCi<ect1,ECCi>=ect2)
        Ci=fct-(0.5*fct/(ect2-ect1))*(ECCi-ect1);
        Ctani=-(0.5*fct/(ect2-ect1));
      elseif and(ECCi<ect2,ECCi>=ect3)
        Ci=.5*fct-(0.5*fct/(ect3-ect2))*(ECCi-ect2);
        Ctani=-(0.5*fct/(ect3-ect2));
      elseif ECCi<ect3
        Ci=0;
        Ctani=0;
     end
     Cci(z)=b*T*h*Ci;% unconfined below tangent concrete force
     Cctani(z)=b*T*h*Ctani;% unconfined below tangent concrete force
     end
    end
       FsTOTAL=sum(Fsi);CcTOTAL=sum(Cci);FssTOTAL=sum(Fssi);Ai=CcTOTAL+FsTOTAL+FssTOTAL-Ptarget;%+FssTOTAL
       FsTOTAL_tan=sum(Fstani);CcTOTAL_tan=sum(Cctani);FssTOTAL_tan=sum(Fstansi);A_tani=CcTOTAL_tan+FsTOTAL_tan+FssTOTAL_tan;%+FssTOTAL_tan
        dxi = A_tani^-1 *(-Ai);
        residual = abs(dxi); % evaluate residual
        it = it + 1; % increment iteration count
        ECCi = ECCi+dxi; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',it,ECCi,residual)
          disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
if it < itermax% iteration control
        fprintf('(+)It is converged in %1.0f iterations - Initial axial strain: %1.6f\n',it,ECCi)
end
else;ECCi=0;end
EC=[abs(ect1) abs(ect2) .33*abs(ect3) .67*abs(ect3) .8*abs(ect3) .9*abs(ect3) abs(ect3) .4*ecuc .5*ecuc .6*ecuc .7*ecuc .8*ecuc .9*ecuc ecuc];q=size(EC,2);
% EC=ECCi+0.00001:0.0001:.0029;q=size(EC,2);
for j=1:q;
    ECC=EC(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
    eC=((ECC*x)/(x-cover));% top Unconfined compression strain
             for u=1:An;
             es=ECCi+(eC*(x-d(u))/x);
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
            if stp==1;
               fs(u)=fy+Esh*(es-esh);
               fstan(u)=(Esh*eC*d(u))/(x)^2;
            elseif stp==2;
               fs(u)=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
               fstan(u)=(2*eC*d(u)*(fu-fy)*(((eC*d(u))/x)+esu-eC))/(x^2*(esu-esh)^2);
            end
       elseif and(es<=(-esh),es>(-esu))
           if stp==1;
              fs(u)=-fy-Esh*(abs(es)-esh);
              fstan(u)=(Esh*eC*d(u))/(x)^2;
           elseif stp==2;
              fs(u)=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
              fstan(u)=(2*eC*d(u)*(fu-fy)*(((eC*d(u))/x)+esu-eC))/(x^2*(esu-esh)^2);
           end
       elseif or(es>=esu,es<=(-esu))
        fs(u)=0;
        fstan(u)=0;
     end
        Fs(u)=As(u)*fs(u);
        Fstan(u)=As(u)*fstan(u);% tangent steel force
        if  and(es>=0,es<ecuc)% in this step: confined concrete force in rebar area is omitted (F=As*fc)
        C_cons(u)=(fcc*r*(es/ecc))/(r-1+(es/ecc)^r);
         c1s=-(d(u)*fcc*r^2*(x-d(u))*es^2*(es*((x-d(u))/(ecc*x)))^(r-1))/(ecc^2*x^3*(((es*(x-d(u)))/(ecc*x))^r +r-1)^2);
         c2s=-(es*fcc*r*(x-d(u)))/(ecc*x^2*(((es*(x-d(u)))/(ecc*x))^r +r-1));
         c3s=(es*fcc*r)/(ecc*x*(((es*(x-d(u)))/(ecc*x))^r +r-1));
         Ctan_cons(u)=c1s+c2s+c3s;
      elseif es>ecuc
         C_cons(u)=0;
         Ctan_cons(u)=0;
       elseif and(es<0,es>=ect1)
        C_cons(u)=-(0.5*Ec*es);
        Ctan_cons(u)=-((0.5*Ec*d(u)*eC)/x^2);
      elseif and(es<ect1,es>=ect2)
        C_cons(u)=+(fct-(0.5*fct/(ect2-ect1))*(es-ect1));
        Ctan_cons(u)=+(0.5*fct*eC*d(u))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        C_cons(u)=-(.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2));
        Ctan_cons(u)=+(0.5*fct*eC*d(u))/((ect3-ect2)*x^2);
      elseif es<ect3
        C_cons(u)=0;
        Ctan_cons(u)=0;  
        end
        Fss(u)=-As(u)*C_cons(u);
        Fstans(u)=-As(u)*Ctan_cons(u);% tangent Minus of concrete force
             end
    for z=1:N;% in this step: concrete force for each fiber is calculated
     %-------------- Cc --------------%
     ec=ECCi+(eC*(x-c(z))/x);
     if and(z>=1,z<=roundn(N*(cover/h),0)) % TOP SECTIONS
     if and(ec>=0,ec<ecu) % unconfined
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
     elseif and(z>roundn((cover/h)*N,0),z<=roundn(N*(1-(cover/h)),0)) % MIDDLE SECTIONS
         if and(ec>=0,ec<ecu) % unconfined
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
      if  and(ec>=0,ec<ecuc)% confined
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
     elseif and(z>roundn((1-(cover/h))*N,0),z<=N) % BELOW SECTIONS
         if and(ec>=0,ec<ecu) % unconfined
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
    for u=1:An;e(u)=x-d(u);end % distance of each rebar from neutral axis
    Pc1=x-.5*h;% Distance from applied force from neutral axis
     for k=1:N;cc(k)=x-c(k);end %distance of each concrete fiber from neutral axis
       FsTOTAL=sum(Fs);CcTOTAL=sum(Cc);FssTOTAL=sum(Fss);A=CcTOTAL+FsTOTAL+FssTOTAL-Ptarget;%+FssTOTAL
       FsTOTAL_tan=sum(Fstan);CcTOTAL_tan=sum(Cctan);FssTOTAL_tan=sum(Fstans);A_tan=CcTOTAL_tan+FsTOTAL_tan+FssTOTAL_tan;%+FssTOTAL_tan
        dx = real(A_tan^-1) *real(-A);
        residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eC,residual)
          disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
       %if abs(x)>h;x=0.5*h;elseif x<0;x=0.5*h;end
    end
    %if it == itermax;break;end % stop the analysis at all because last Convergence     
     
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.6f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eC,x,(eC/x)*1000,(Fs*e'+Cc*cc'+Pc1*Ptarget)*.000001)
        end
TCCS(j)=(fcc*r*(ECC/ecc))/(r-1+(ECC/ecc)^r);% Top Confined compression stress
if and(eC>0,eC<ecu);ECCC(j)=eC;TUCS(j)=(fc*R*(eC/ec0))/(R-1+(eC/ec0)^R);elseif and(eC>=ecu,eC<esp);ECCC(j)=eC;TUCS(j)=fcu*(1-((eC-ecu)/(esp-ecu)));elseif eC>=esp;ECCC(j)=eC;TUCS(j)=0;end% Top Unconfined compression stress
TSCSs(j)=eC*(x-d(1))/x;TSCS(j)=fs(1);% Top Steel compression strain-stress
BSCSs(j)=eC*(x-d(6))/x;BSCS(j)=fs(6);% Bottom Steel compression strain-stress
if CD==0;CrackDepth(j)=0; else CrackDepth(j)=max(abs(CD));end;% Crack Depth of each increment
% Calculate Moment and Curavture
Cur(j)=(eC/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=roundn((Fs*e'+Fss*e'+Cc*cc'-Ptarget*Pc1),-3)*.000001;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
Cur=[0 Cur];Mom=[0 Mom];
s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));end % Flextural Rigidity
if eC == ecuc;fprintf('\n       ## Confined Concrete Strain Reached to Ultimate Strain: %1.4f ## \n\n',eC);end
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
0.002542
0.005083
0.007625
0.01017
0.01271
0.01525
0.01779
0.02033
0.02287
0.02542
0.02829
0.03116
0.03404
0.03691
0.03978
0.04265
0.04553
0.0484
0.05127
0.05415
0.05702
0.05989
0.06277
0.06564
0.06851
0.07139
0.07426
0.07713
0.08001
0.08288
0.08288];

MomXT=[0
137.4
212
277.6
337.2
391.7
395.7
397.9
399.4
400.5
401.1
401.2
401.1
401
400.3
399.2
395.7
381.1
367
360.4
356.2
355.3
355
354.5
353.7
353.8
353.6
353.3
352.8
352.3
351.8
351.8];

%% imaging
figure(1)
IMAGE1=imread('FiberConfinedConcreteMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE2=imread('FiberConfinedConcreteMomentCurvature-image2.jpg');
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
title(' Plotting of the Top concrete strain and Crack depth ','Color','b');grid on;
figure(6)
plot(EC,XX,'--','LineWidth',2);xlabel('Top confined concrete fiber strain');ylabel('Neutral axis (mm)')
title(' Plotting of the Top confined concrete strain and neutral axis ','Color','b');grid on
figure(7)
plot(EC,CUR,'--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top concrete strain and Curvature ','Color','b');grid on
figure(8)
plot(EC,TCCS,ECCC,TUCS,'g--','LineWidth',3);xlabel('Top fiber concrete strain (\epsilon_c)');ylabel('Top fiber concrete stress (N/mm^2)')
title(' Plotting of the Top fiber confined and unconfined concrete strain-stress during the analysis','Color','b');
legend('Confined','Unconfined','Location','NorthEastOutside');grid on
figure(9)
plot(TSCSs,TSCS,BSCSs,BSCS,'g--','LineWidth',3);xlabel('Steel strain (\epsilon_s)');ylabel('Steel stress (N/mm^2)')
title(' Plotting of the Top  and Bottom steel strain-stress during the analysis','Color','b');
legend('Top steel','Bottom steel','Location','NorthEastOutside');grid on
figure(10)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(11)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(12)
plot(CUR,EI,'black--','LineWidth',3)
title('# Flextural Rigidity(EI)-Curvature diagram #','Color','b');
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(13)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral axis-Flextural Rigidity(EI) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(14)
plot(Cur,Mom,X,Y,'r--',CurXT,MomXT,'g-.','LineWidth',3)
title(['# CONFINED SECTION MOMENT CURVATURE ANALYSIS DIAGRAM #','     Target axial load : ',int2str(Ptarget*.001),' (kN)  -  EI : ',int2str(Y(2)/X(2)),' (kN.m^2)  -  EI_t : ',int2str((Y(3)-Y(2))/(X(3)-X(2))),' (kN.m^2)  -  Ductility Rito : ',int2str(Ductility_Rito)],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Confined','Confined bilinear fitted','XTRACT','Location','NorthEastOutside');grid on
%% Output Data to .txt file
fid = fopen('FiberConfinedConcreteMomentCurvatureAxialLoad-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*   Moment-Curvature Analysis of Confined concrete section  *\n');
fprintf(fid,'*   With Axial Load effect                                  *\n');
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
fprintf(fid,'+===============================================================================+\n');
fclose(fid);