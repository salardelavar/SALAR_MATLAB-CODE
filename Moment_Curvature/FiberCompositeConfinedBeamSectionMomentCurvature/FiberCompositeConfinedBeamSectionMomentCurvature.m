%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
%  Moment-Curvature Analysis of Composite Confined Beam     %
%  Section                                                  %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given:Section Properties , SteelSection properties ,       %
%Calculate: Moment-Curavture                                %
% Newton-Raphson Method : Tangent procedure                 %
%***********************************************************%
%   _    ______________________________________             %
%   |   |                                      |            %
%       |     #     #     #     #    #    #    |            %
%       |     #                           #    |            %
% Beff  |    As1   As2   As3   As4  As5  As6   |            %
%       |     #                           #    |            %
%   |   |     #     #     #     #    #    #    |            %
%   _   |______________________________________|            %
%       |<-               Teff               ->|            %
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
%        _     __                 __                        %   
%        ^    |  |               |  |                       %
%        |    |  |               |  |                       %
%             |  |               |  |                       %
%             |  |_______________|  |  +                    %
%       bf    |   _______________   |  tw                   %
%             |  |               |  |  +                    %
%        |    |  |               |  |                       %
%        v    |  |               |  |                       %
%        -     ---               ---                        %
%             |<-   tf1+tf2+ hw   ->|                       % 
%    X                                                      %
%    ^                                                      %
%    |             (Moment - Curvature along X axis)        %
%    |                                                      %
%    +----> Y                                               %
%***********************************************************%
clear all;close all;clc;
%% Section Properties
Beff=1000; % Width of deck
Teff=100; % Depth of deck
%As:  As1      As2
As=[1000 1000];
%d:d1  d2 
d=[30 60];
cover=10; % Concrete cover - NOTE: Cover must less than d(1)
tf1=9.2;% [mm] I section thickness on Top flange
bf1=110;% [mm] I section width on Top flange
tw=5.9;% [mm] I section thickness of Web
hw=201.6;% [mm] Height of web
tf2=9.2;% [mm] I section thickness on Bottom flange
bf2=110;% [mm] I section width on Bottom flange
ptf2=10;% [mm] Plate section thickness on Bottom flange
pbf2=80;% [mm] Plate section width on Bottom flange
h=Teff+tf1+tf2+hw+ptf2;% [mm] Height of Section
N=h*2;% Number of steel section Fiber 
itermax = 1000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
x=.25*h;% initial guess of Neuteral axis
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
esp=ecu+.002; % Spalling concrete strain
Ec=5000*sqrt(fc);
ec0=(2*fc)/Ec;
fct=-0.7*sqrt(fc);% Concrete tension stress
ect1=(2*fct)/Ec;ect2=(2.625*fct)/Ec;ect3=(9.292*fct)/Ec;% Concrete tension strain
%% Steel Reinforcing Properties
fys =400;% Yield strength of steel reinforcing (N/mm^2)
Ess =2e5;% Modulus of elasticity of steel (N/mm^2)
fus=1.5*fys;% Ultimate steel stress
eys=fys/Ess;% Yeild steel strain
eshs=0.01;% Strain at steel strain-hardening
esus=0.09;% Ultimate steel strain
stp= 2;% Kind of strain hardening:  1:linear   2:curve
Eshs=(fus-fys)/(esus-eshs);
%% Reinforcing stirrup Properties
fyh =300;% [MPa] Transverse Reinforcing Bar (stirrup) Yield Stress
Nx=2;% Total of transverse hoop legs in the X
Ny=4;% Total of transverse hoop legs in the Y
diastirrup=8;% [mm] Cross-sectional hoop diameter
wi=6.4e+9;% [mm^2] Define total summation of confined reinforcing steel in power.2
s =100;% [mm] Tie Spacing Along Member
%%% monitor cpu time
starttime = cputime;
%% Calculating concrete properties
Asp=3.1415*(diastirrup^2)/4;
cx=Teff-2*cover-diastirrup;%confined core dimensions in the X 
cy=Beff-2*cover-diastirrup;%confined core dimensions in the Y
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
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure ------------------------%
N1=N/10;% Number of top cover concrete section
N2=N/10;% Number of middle concrete section
N3=N/10;% Number of below cover concrete section
N4=N/10;% Number of steel section top flange Fiber  
N5=4*N/10;% Number of steel section  web Fiber
N6=N/10;% Number of steel section bottom flange Fiber
N7=N/10;% Number of steel section bottom Plate flange Fiber 
An=size(As,2);R1=(1/N1);R2=(1/N2);R3=(1/N3);R4=(1/N4);R5=(1/N5);R6=(1/N6);R7=(1/N7);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*cover;end
for k=1:N2;c2(k)=cover+(.5*R2+R2*(k-1))*(Teff-2*cover);end 
for k=1:N3;c3(k)=cover+(Teff-2*cover)+(.5*R3+R3*(k-1))*cover;end 
for k=1:N4;c4(k)=Teff+(.5*R4+R4*(k-1))*tf1;end % distance of each top flange Fiber 
for k=1:N5;c5(k)=Teff+tf1+(.5*R5+R5*(k-1))*hw;end % distance of each web Fiber
for k=1:N6;c6(k)=Teff+tf1+hw+(.5*R6+R6*(k-1))*tf2;end % distance of each bottom flange Fiber
for k=1:N7;c7(k)=Teff+tf1+hw+tf2+(.5*R7+R7*(k-1))*ptf2;end % distance of each bottom flange Fiber
c=[c1 c2 c3 c4 c5 c6 c7];
ES=[.1*abs(ect1) .33*abs(ect1) .67*abs(ect1) abs(ect1) abs(ect2) ...
    .33*abs(ect3) .5*abs(ect3) .67*abs(ect3) .8*abs(ect3) .9*abs(ect3) abs(ect3) ...
    .015*ecuc .02*ecuc .03*ecuc .04*ecuc .05*ecuc .06*ecuc .07*ecuc .08*ecuc .09*ecuc ...
    .1*ecuc];q=size(ES,2); %%.15*ecuc .2*ecuc .3*ecuc .4*ecuc .5*ecuc .6*ecuc .7*ecuc .8*ecuc .9*ecuc ecuc];
%ES=.5*abs(ect1):((.5*abs(ect1)+10^-3)/ecuc)*.5*abs(ect1)+10^-3:0.01833;q=size(ES,2);
x=.5*h;% initial guess of Neuteral axis
for j=1:q;
    eS=ES(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        %---------------- As -------------%
        for u=1:An;
        ess=eS*(x-d(u))/x;
       if and(ess>=0,ess<eys)
        fsa(u)=Ess*ess;
        fstana(u)=(Ess*eS*d(u))/(x)^2;
       elseif and(ess<0,ess>(-eys))
        fsa(u)=Ess*ess;
        fstana(u)=(Ess*eS*d(u))/(x)^2;
       elseif  and(ess>=eys,ess<eshs)
        fsa(u)=fys;
        fstana(u)=0;
       elseif and(ess<=(-eys),ess>(-eshs))
        fsa(u)=-fys;
        fstana(u)=0;
        elseif  and(ess>=eshs,ess<esus)
        fsa(u)=fus-(fus-fys)*((esus-abs(ess))/(esus-eshs))^2;
        fstana(u)=(2*eS*d(u)*(fus-fys)*(((eS*d(u))/x)+esus-eS))/(x^2*(esus-eshs)^2);
       elseif and(ess<=(-eshs),ess>(-esus))
        fsa(u)=-fus+(fus-fys)*((esus-abs(ess))/(esus-eshs))^2;
        fstana(u)=(2*eS*d(u)*(fus-fys)*(((eS*d(u))/x)+esus-eS))/(x^2*(esus-eshs)^2);
       elseif or(ess>=esus,ess<=(-esus))
        fsa(u)=0;
        fstana(u)=0;
       end
        FSA(j,u)=fsa(u);
        FsA(u)=As(u)*fsa(u);
        FstanA(u)=As(u)*fstana(u);% tangent steel force
      if and(ess>=0,ess<=ecuc)% in this step: Confined concrete force in rebar area is omitted (F=As*fc)
         Cs=(fcc*r*(ess/ecc))/(r-1+(ess/ecc)^r);
         cap1=-(d(u)*fcc*r^2*(x-d(u))*ess^2*(ess*((x-d(u))/(ecc*x)))^(r-1))/(ecc^2*x^3*(((ess*(x-d(u)))/(ecc*x))^r +r-1)^2);
         cap2=-(ess*fcc*r*(x-d(u)))/(ecc*x^2*(((ess*(x-d(u)))/(ecc*x))^r +r-1));
         cap3=(ess*fcc*r)/(ecc*x*(((ess*(x-d(u)))/(ecc*x))^r +r-1));
         Ctans=cap1+cap2+cap3;
      elseif ess > ecuc
         Cs=0; 
         Ctans=0;
      elseif and(ess<0,ess>=ect1)
        Cs=-0.5*Ec*ess;
        Ctans=-(0.5*Ec*d(u)*eS)/x^2;
      elseif and(ess<ect1,ess>=ect2)
        Cs=fct+(0.5*fct/(ect2-ect1))*(ess-ect1);
        Ctans=+(0.5*fct*eS*d(u))/((ect2-ect1)*x^2);
      elseif and(ess<ect2,ess>=ect3)
        Cs=-(.5*fct-(0.5*fct/(ect3-ect2))*(ess-ect2));
        Ctans=+(0.5*fct*eS*d(u))/((ect3-ect2)*x^2);
      elseif ess<ect3
        Cs=0;
        Ctans=0;
        end
        CS(u)=Cs;CtanS(u)=Ctans;
        Fss(u)=-As(u)*CS(u);
        Fstans(u)=-As(u)*CtanS(u);% tangent Minus of concrete force
        end
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c(z))/x;            
        if and(z>=1,z<=roundn(N/10,0)) % TOP UNCONFINED CONCRETE SECTION
      if and(es>=0,es<=ecu) % unconfined
         fs=(fc*R*(es/ec0))/(R-1+(es/ec0)^R);
         Ca1=-(c(z)*fc*R^2*(x-c(z))*es^2*(es*((x-c(z))/(ec0*x)))^(R-1))/(ec0^2*x^3*(((es*(x-c(z)))/(ec0*x))^R +R-1)^2);
         Ca2=-(es*fc*R*(x-c(z)))/(ec0*x^2*(((es*(x-c(z)))/(ec0*x))^R +R-1));
         Ca3=(es*fc*R)/(ec0*x*(((es*(x-c(z)))/(ec0*x))^R +R-1));
         fstan=Ca1+Ca2+Ca3;
      elseif and(es>ecu,es<esp)
         fs=fcu*(1-((es-ecu)/(esp-ecu)));
         fstan=-(eS*c(z)*fcu)/((esp-ecu)*x^2);
      elseif es >= esp
         fs=0; 
         fstan=0;
      elseif and(es<0,es>=ect1)
        fs=0.5*Ec*es;
        fstan=(0.5*Ec*c(z)*eS)/x^2;
      elseif and(es<ect1,es>=ect2)
        fs=fct-(0.5*fct/(ect2-ect1))*(es-ect1);
        fstan=-(0.5*fct*eS*c(z))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        fs=.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2);
        fstan=-(0.5*fct*eS*c(z))/((ect3-ect2)*x^2);
      elseif es<ect3
        fs=0;
        fstan=0;
     end
      Fs(z)=Beff*R1*cover*fs;
      Fstan(z)=Beff*R1*cover*fstan;% tangent concrete force
      elseif and(z>=roundn(N/10,0),z<=roundn(2*N/10,0)) % MIDDLE CONFINED CONCRETE SECTION
      if  and(es>=0,es<=ecuc)% confined
        fs=(fcc*r*(es/ecc))/(r-1+(es/ecc)^r);
         ca1=-(c(z)*fcc*r^2*(x-c(z))*es^2*(es*((x-c(z))/(ecc*x)))^(r-1))/(ecc^2*x^3*(((es*(x-c(z)))/(ecc*x))^r +r-1)^2);
         ca2=-(es*fcc*r*(x-c(z)))/(ecc*x^2*(((es*(x-c(z)))/(ecc*x))^r +r-1));
         ca3=(es*fcc*r)/(ecc*x*(((es*(x-c(z)))/(ecc*x))^r +r-1));
         fstan=ca1+ca2+ca3;
      elseif es>ecuc
         fs=0;
         fstan=0;
       elseif and(es<0,es>=ect1)
        fs=0.5*Ec*es;
        fstan=(0.5*Ec*c(z)*eS)/x^2;
      elseif and(es<ect1,es>=ect2)
        fs=fct-(0.5*fct/(ect2-ect1))*(es-ect1);
        fstan=-(0.5*fct*eS*c(z))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        fs=.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2);
        fstan=-(0.5*fct*eS*c(z))/((ect3-ect2)*x^2);
      elseif es<ect3
        fs=0;
        fstan=0;  
      end
     Fs(z)=Beff*R2*(Teff-2*cover)*fs; %confined and unconfined middle concrete force
     Fstan(z)=Beff*R2*(Teff-2*cover)*fstan;%confined and unconfined middle tangent concrete force
        elseif and(z>=roundn(2*N/10,0),z<=roundn(3*N/10,0)) % TOP UNCONFINED CONCRETE SECTION
      if and(es>=0,es<=ecu) % unconfined
         fs=(fc*R*(es/ec0))/(R-1+(es/ec0)^R);
         Ca1=-(c(z)*fc*R^2*(x-c(z))*es^2*(es*((x-c(z))/(ec0*x)))^(R-1))/(ec0^2*x^3*(((es*(x-c(z)))/(ec0*x))^R +R-1)^2);
         Ca2=-(es*fc*R*(x-c(z)))/(ec0*x^2*(((es*(x-c(z)))/(ec0*x))^R +R-1));
         Ca3=(es*fc*R)/(ec0*x*(((es*(x-c(z)))/(ec0*x))^R +R-1));
         fstan=Ca1+Ca2+Ca3;
      elseif and(es>ecu,es<esp)
         fs=fcu*(1-((es-ecu)/(esp-ecu)));
         fstan=-(eS*c(z)*fcu)/((esp-ecu)*x^2);
      elseif es >= esp
         fs=0; 
         fstan=0;
      elseif and(es<0,es>=ect1)
        fs=0.5*Ec*es;
        fstan=(0.5*Ec*c(z)*eS)/x^2;
      elseif and(es<ect1,es>=ect2)
        fs=fct-(0.5*fct/(ect2-ect1))*(es-ect1);
        fstan=-(0.5*fct*eS*c(z))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        fs=.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2);
        fstan=-(0.5*fct*eS*c(z))/((ect3-ect2)*x^2);
      elseif es<ect3
        fs=0;
        fstan=0;
     end
      Fs(z)=Beff*R3*cover*fs;
      Fstan(z)=Beff*R3*cover*fstan;% tangent concrete force
        elseif and(z>=roundn(3*N/10,0),z<=roundn(4*N/10,0)) % TOP PLATE FLANGE    
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
     Fs(z)=bf1*R4*tf1*fs;% force on top flange
     Fstan(z)=bf1*R4*tf1*fstan;% tangent steel force
        elseif and(z>roundn(4*N/10,0),z<=roundn(8*N/10,0)) % MIDDLE WEB
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
     Fs(z)=tw*R5*hw*fs;% force on web
     Fstan(z)=tw*R5*hw*fstan;% tangent steel force
        elseif and(z>roundn(8*N/10,0),z<=roundn(9*N/10,0)) % BOTTOM I FLANGE
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
      Fs(z)=bf2*R6*tf2*fs;% force on bottom flange
      Fstan(z)=bf2*R6*tf2*fstan;% tangent steel force
     elseif and(z>roundn(9*N/10,0),z<=N) % BOTTOM PLATE FLANGE
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
      Fs(z)=pbf2*R7*ptf2*fs;% force on plate bottom flange
      Fstan(z)=pbf2*R7*ptf2*fstan;% tangent steel force
        end 
     Ss(z)=es;SS(j,z)=Ss(z);% sction Fiber Strain
     CFS(j,z)=fs;% sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);FsATOTAL=sum(FsA);FssTOTAL=sum(Fss);A=FsTOTAL+FsATOTAL+FssTOTAL;
    FsTOTAL_tan=sum(Fstan);FsATOTAL_tan=sum(FstanA);FssTOTAL_tan=sum(Fstans);A_tan=FsTOTAL_tan+FsATOTAL_tan+FssTOTAL_tan;
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
    for u=1:An;e(u)=x-d(u);end % distance of each rebar from neutral axis        
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc(k)=x-c(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc'+FsA*e'+Fss*e')*10^-6)
        end                        
TSCSs(j)=eS*(x-d(1))/x;TSCS(j)=fsa(1);% Top Steel compression strain-stress
BSCSs(j)=eS*(x-d(2))/x;BSCS(j)=fsa(An);% Bottom Steel compression strain-stress       
eST=eS*(x-Teff)/x;if and(eST>0,eST<ey);EST(j)=eST;TUST(j)=Es*eST;elseif and(eST>=ey,eST<esh);EST(j)=eST;TUST(j)=fy;elseif and(eST>=esh,eST<esu);EST(j)=eST;TUST(j)=fu-(fu-fy)*((esu-abs(eST))/(esu-esh))^2;elseif and(eST<0,eST>(-ey));EST(j)=eST;TUST(j)=Es*eST;elseif and(eST<=(-ey),eST>(-esh));EST(j)=eST;TUST(j)=-fy;elseif and(eST<=(-esh),eST>(-esu));EST(j)=eST;TUST(j)=-fu+(fu-fy)*((esu-abs(eST))/(esu-esh))^2;end% Top Steel section compression stress
eSB=eS*(x-h)/x;if and(eSB>0,eSB<ey);ESB(j)=eSB;TUSB(j)=Es*eSB;elseif and(eSB>=ey,eSB<esh);ESB(j)=eSB;TUSB(j)=fy;elseif and(eSB>=esh,eSB<esu);ESB(j)=eSB;TUSB(j)=fu-(fu-fy)*((esu-abs(eSB))/(esu-esh))^2;elseif and(eSB<0,eSB>(-ey));ESB(j)=eSB;TUSB(j)=Es*eSB;elseif and(eSB<=(-ey),eSB>(-esh));ESB(j)=eSB;TUSB(j)=-fy;elseif and(eSB<=(-esh),eSB>(-esu));ESB(j)=eSB;TUSB(j)=-fu+(fu-fy)*((esu-abs(eSB))/(esu-esh))^2;end% Bottom Steel section compression stress
BCTS(j)=(eS*(h-cover-x))/(x-cover);% Bottom Confined tension strain
TCCS(j)=(fcc*r*(eS/ecc))/(r-1+(eS/ecc)^r);% Top Confined compression stress
if and(eS>0,eS<ecu);ECCC(j)=eS;TUCS(j)=(fc*R*(eS/ec0))/(R-1+(eS/ec0)^R);elseif and(eS>=ecu,eS<esp);ECCC(j)=eS;TUCS(j)=fcu*(1-((eS-ecu)/(esp-ecu)));elseif eS>=esp;ECCC(j)=eS;TUCS(j)=0;end% Top Unconfined compression stress 
% Calculate Moment and Curavture
Cur(j)=(eS/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*cc'+FsA*e'+Fss*e')*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
Cur=[0 Cur];Mom=[0 Mom];
%s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));end % Flextural Rigidity
s=size(Cur,2);for i=1:s-1;EI(i)=Mom(i)/Cur(i);end % Flextural Rigidity
if eS == ecuc;fprintf('\n      ##  Concrete Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
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
disp('=    Section curve fitted     =');
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
%% imaging
figure(1)
IMAGE1=imread('FiberCompositeConfinedBeamSectionMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE1=imread('FiberCompositeConfinedBeamSectionMomentCurvature-image2.jpg');
image(IMAGE1);axis image;axis off;
%% Plot
figure(3)
plot(TSCSs,TSCS,BSCSs,BSCS,'g--','LineWidth',3);xlabel('Steel strain (\epsilon_s)');ylabel('Steel stress (N/mm^2)')
title(' Plotting of the Top  and Bottom steel rebar strain-stress during the analysis','Color','b');
legend('Top steel rebar','Bottom steel rebar','Location','NorthEastOutside');grid on;
figure(4)
P1=plot(EST,TUST,ESB,TUSB,'r--');set(P1,'LineWidth',3);
xlabel('Strain');ylabel('Stress (kN/mm^2)')
legend('Top steel fiber','Bottom steel fiber','Location','NorthEastOutside');grid on;
title('Top and Bottom Steel Fiber Section Stress-Strain Diagram','color','b')
figure(5)
plot(ES,TCCS,ECCC,TUCS,'g--','LineWidth',3);xlabel('Top fiber concrete strain ');ylabel('Top fiber concrete stress (N/mm^2)')
title(' Plotting of the Top fiber confined and unconfined concrete strain-stress during the analysis','Color','b');
legend('Confined','Unconfined','Location','NorthEastOutside');grid on
figure(6)
p1=plot(I,DX,'b--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
figure(7)
p1=plot(I,IT,'b--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
figure(8)
plot(ES,XX,'--','LineWidth',3);xlabel('Top Concrete strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top concrete strain and Neuteral axis ','Color','b');grid on
figure(9)
plot(ES,CUR,'--','LineWidth',3);xlabel('Top concrete strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top concrete strain and Curvature ','Color','b');grid on
figure(10)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the steel and concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(11)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the steel and concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
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
plot(Cur,Mom,X,Y,'r--','LineWidth',3)
title(['# COMPOSITE BEAM CONFINED SECTION MOMENT-CURVATURE DIAGRAM #','  EI_e_l_a : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Composite','Composite - bilinear fitted','Location','NorthEastOutside');grid on;
figure(15)
Cur02=[0
    0.0003
    0.0010
    0.0021
    0.0031
    0.0013
    0.0027
    0.0040
    0.0047
    0.0118
    0.0151
    0.0177
    0.0204
    0.0183
    0.0269
    0.0372
    0.0475
    0.0572
    0.0667
    0.0760
    0.0850];
Mom02=[0
    9.4411
   30.9767
   62.3519
   92.2741
   40.5514
   81.3976
  120.1220
  139.4997
  184.7034
  190.9853
  194.9387
  198.4406
  195.7591
  205.0244
  207.2787
  208.5784
  209.5055
  210.2265
  210.8218
  211.3304];
plot(Cur,Mom,Cur02,Mom02,'--r','LineWidth',3)
title(['# COMPOSITE BEAM CONFINED AND UNCONFINED SECTION MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Confined Composite Section','Unconfined Composite Section','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberCompositeConfinedBeamSectionMomentCurvature-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'* Moment-Curvature Analysis of Composite Confined Beam      *\n');
fprintf(fid,'* Section                                                   *\n');
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
fprintf(fid,'*       |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'* Beff  |    As1   As2   As3   As4  As5  As6   |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*   |   |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*   _   |______________________________________|            *\n');
fprintf(fid,'*       |<-               Teff               ->|            *\n');
fprintf(fid,'*       |<-d1->|                                            *\n');
fprintf(fid,'*       |<-  d2   ->|                                       *\n');
fprintf(fid,'*       |<-     d3      ->|                                 *\n');
fprintf(fid,'*       |<-        d4          ->|                          *\n');
fprintf(fid,'*       |<-            d5          ->|                      *\n');
fprintf(fid,'*       |<-               d6             >|                 *\n');
fprintf(fid,'*    X                                                      *\n');
fprintf(fid,'*    ^                                                      *\n');
fprintf(fid,'*    |             (Moment - Curvature along X axis)        *\n');
fprintf(fid,'*    |                                                      *\n');
fprintf(fid,'*    +----> Y                                               *\n');
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
fprintf(fid,'\n');
fprintf(fid,'+==============================+\n');
fprintf(fid,'=   Section Moment-Curvature   =\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur;Mom]);
fprintf(fid,'+=============================+\n');
fprintf(fid,'=    Section curve fitted     =\n');
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