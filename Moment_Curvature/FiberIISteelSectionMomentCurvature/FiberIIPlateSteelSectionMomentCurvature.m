%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
%   Moment-Curvature Analysis of Double I steel sections    %
%   with Double Plates on Flange                            %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given: Steel section properties                            %
%Calculate: Moment-Curavture                                %
% Newton-Raphson Method : Tangent procedure                 %
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
tf1=9.2;% [mm] I section thickness on Top flange
bf1=110;% [mm] I section width on Top flange
tw=5.9;% [mm] I section thickness of Web
hw=201.6;% [mm] Height of web
tf2=9.2;% [mm] I section thickness on Bottom flange
bf2=110;% [mm] I section width on Bottom flange
ptf1=10;% [mm] Plate section thickness on Top flange
pbf1=150;% [mm] Plate section width on Top flange
ptf2=10;% [mm] Plate section thickness on Bottom flange
pbf2=150;% [mm] Plate section width on Bottom flange
h=ptf1+ptf2+tf1+tf2+hw;% [mm] Height of Section
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
tolerance = 1e-9; % specified tolerance for convergence
x=.25*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure ------------------------%
N1=N/10;% Number of I steel section top plate flange Fiber
N2=N/10;% Number of steel section top flange Fiber 
N3=6*N/10;% Number of I steel section web Fiber
N4=N/10;% Number of I steel section bottom flange Fiber
N5=N/10;% Number of steel section bottom Plate flange Fiber 
R1=(1/N1);R2=(1/N2);R3=(1/N3);R4=(1/N4);R5=(1/N5);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*ptf1;end % distance of each top Plate flange Fiber 
for k=1:N2;c2(k)=ptf1+(.5*R2+R2*(k-1))*tf1;end % distance of each top I flange Fiber
for k=1:N3;c3(k)=ptf1+tf1+(.5*R3+R3*(k-1))*hw;end % distance of each web Fiber
for k=1:N4;c4(k)=ptf1+tf1+hw+(.5*R4+R4*(k-1))*tf2;end % distance of each bottom I flange Fiber
for k=1:N5;c5(k)=ptf1+tf1+hw+tf2+(.5*R5+R5*(k-1))*ptf2;end % distance of each bottom Plate flange Fiber
c=[c1 c2 c3 c4 c5];
ES=ey:.001:esu;q=size(ES,2);
for j=1:q;
    eS=ES(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c(z))/x;            
        if and(z>=1,z<=roundn(N/10,0)) % Top Plate of flange
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
     Fs(z)=pbf1*R1*ptf1*fs;% force on top flange
     Fstan(z)=pbf1*R1*ptf1*fstan;% tangent steel force
        elseif and(z>roundn(N/10,0),z<=roundn(2*N/10,0)) % Top flange of I
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
     Fs(z)=2*bf1*R2*tf1*fs;% force on top flange
     Fstan(z)=2*bf1*R2*tf1*fstan;% tangent steel force
        elseif and(z>roundn(2*N/10,0),z<=roundn(8*N/10,0)) % Middle web
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
     Fs(z)=2*tw*R3*hw*fs;% force on web
     Fstan(z)=2*tw*R3*hw*fstan;% tangent steel force
     elseif and(z>roundn(8*N/10,0),z<=roundn(9*N/10,0)) % Below flange of I
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
      Fs(z)=2*bf2*R4*tf2*fs;% force on bottom flange
      Fstan(z)=2*bf2*R4*tf2*fstan;% tangent steel force
      elseif and(z>roundn(9*N/10,0),z<=N) % below Plate of flange
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
      Fs(z)=pbf2*R5*ptf2*fs;% force on bottom flange
      Fstan(z)=pbf2*R5*ptf2*fstan;% tangent steel force
        end
     Ss(z)=es;SS(j,z)=Ss(z);%steel sction Fiber Strain
     CFS(j,z)=fs;%steel sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);A=FsTOTAL;
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
     for k=1:N;cc(k)=x-c(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc')*10^-6)
        end
 if and(eS>0,eS<ey);ESS(j)=eS;TUSS(j)=Es*eS;elseif and(eS>=ey,eS<esh);ESS(j)=eS;TUSS(j)=fy;elseif and(eS>=esh,eS<esu);ESS(j)=eS;TUSS(j)=fu-(fu-fy)*((esu-abs(eS))/(esu-esh))^2;end% Top Steel section compression stress        
% Calculate Moment and Curavture
Cur(j)=(eS/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*cc')*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
if roundn(eS,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
Cur=[0 Cur];Mom=[0 Mom];
s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));end % Flextural Rigidity
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
%% Elastic-Perfect Plastic Steel Section Moment-Curvature Analysis (simple)
%  in this simple analysis : bf1=bf2 & tf1=tf2 & tw1=tw2 & hw1=hw2 & pbf1=pbf2 & ptf1=ptf2
Ie1=2*(((tw*hw^3)/12)+2*(((bf1*tf1^3)/12)+(bf1*tf1)*(0.5*h-0.5*tf1-ptf1)^2)+(((pbf1*ptf1^3)/12)+(pbf1*ptf1)*(0.5*h-0.5*ptf1)^2));
phiY1=ey/(.5*h);phiU1=esu/(.5*h);
My1=(fy*Ie1)/(.5*h); % Yeild Moment
Mu1=fy*[2*(2*bf1*tf1*(.5*h-.5*tf1)+pbf1*ptf1*(.5*h-.5*ptf1)+(tw*hw^2)/4)];  % Ultimate Moment
EIe1=(My1/phiY1)*1e-9;% Elastic Flextural Rigidity
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
Material_Ductility_Rito=esu/esh;
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
%% C++ outout FiberRecSectionMomentCurvatureLayerStrainSteel
Curvature=[0
7.29167e-006
1.45833e-005
2.1875e-005
2.91667e-005
3.64583e-005
4.375e-005
5.10417e-005
5.83333e-005
6.5625e-005
7.29167e-005
8.02083e-005
8.75e-005
9.47917e-005
0.000102083
0.000109375
0.000116667
0.000123958
0.00013125
0.000138542
0.000145833
0.000153125
0.000160417
0.000167708
0.000175
0.000182292
0.000189583
0.000196875
0.000204167
0.000211458
0.00021875
0.000226042
0.000233333
0.000240625
0.000247917
0.000255208
0.0002625
0.000269792
0.000277083
0.000284375
0.000291667
0.000298958
0.00030625
0.000313542
0.000320833
0.000328125
0.000335417
0.000342708
0.00035
0.000357292
0.000364583
0.000371875
0.000379167
0.000386458
0.00039375
0.000401042
0.000408333
0.000415625
0.000422917
0.000430208
0.0004375
0.000444792
0.000452083
0.000459375
0.000466667
0.000473958
0.00048125
0.000488542
0.000495833
0.000503125
0.000510417
0.000517708
0.000525
0.000532292
0.000539583
0.000546875
0.000554167
0.000561458
0.00056875
0.000576042
0.000583333
0.000590625
0.000597917
0.000605208
0.0006125
0.000619792
0.000627083
0.000634375
0.000641667
0.000648958
0.00065625
0.000663542
0.000670833
0.000678125
0.000685417
0.000692708
0.0007
0.000707292
0.000714583
0.000721875
0.000729167
0.000736458
0.00074375
0.000751042
0.000758333
0.000765625
0.000772917
0.000780208
0.0007875
0.000794792
0.000802083
0.000809375
0.000816667
0.000823958
0.00083125
0.000838542
0.000845833
0.000853125
0.000860417
0.000867708
0.000875
0.000882292
0.000889583
0.000896875
0.000904167
0.000911458
0.00091875
0.000926042
0.000933333
0.000940625
0.000947917
0.000955208
0.0009625
0.000969792
0.000977083
0.000984375
0.000991667
0.000998958
0.00100625
0.00101354
0.00102083
0.00102812
0.00103542
0.00104271
0.00105
0.00105729
0.00106458
0.00107187
0.00107917
0.00108646
0.00109375
0.00110104
0.00110833
0.00111562
0.00112292
0.00113021
0.0011375
0.00114479
0.00115208
0.00115937
0.00116667
0.00117396
0.00118125
0.00118854
0.00119583
0.00120312
0.00121042
0.00121771
0.001225
0.00123229
0.00123958
0.00124687
0.00125417
0.00126146
0.00126875
0.00127604
0.00128333
0.00129062
0.00129792
0.00130521
0.0013125
0.00131979
0.00132708
0.00133437
0.00134167
0.00134896
0.00135625
0.00136354
0.00137083
0.00137812
0.00138542
0.00139271
0.0014
0.00140729
0.00141458
0.00142187
0.00142917
0.00143646
0.00144375
0.00145104
0.00145833
0.00146562
0.00147292
0.00148021
0.0014875
0.00149479
0.00150208
0.00150937
0.00151667
0.00152396
0.00153125
0.00153854
0.00154583
0.00155312
0.00156042
0.00156771
0.001575
0.00158229
0.00158958
0.00159687
0.00160417
0.00161146
0.00161875
0.00162604
0.00163333
0.00164062
0.00164792
0.00165521
0.0016625
0.00166979
0.00167708
0.00168437
0.00169167
0.00169896
0.00170625
0.00171354
0.00172083
0.00172812
0.00173542
0.00174271
0.00175
0.00175729
0.00176458
0.00177187
0.00177917
0.00178646
0.00179375
0.00180104
0.00180833
0.00181562
0.00182292
0.00183021
0.0018375
0.00184479
0.00185208
0.00185937
0.00186667
0.00187396
0.00188125
0.00188854
0.00189583
0.00190312
0.00191042
0.00191771
0.001925
0.00193229
0.00193958
0.00194687
0.00195417
0.00196146
0.00196875
0.00197604
0.00198333
0.00199062
0.00199792
0.00200521
0.0020125
0.00201979
0.00202708
0.00203437
0.00204167
0.00204896
0.00205625
0.00206354
0.00207083
0.00207812
0.00208542
0.00209271
0.0021
0.00210729
0.00211458
0.00212187
0.00212917
0.00213646
0.00214375
0.00215104
0.00215833
0.00216562
0.00217292
0.00218021
0.0021875
0.00219479
0.00220208
0.00220937
0.00221667
0.00222396
0.00223125
0.00223854
0.00224583
0.00225312
0.00226042
0.00226771
0.002275
0.00228229
0.00228958
0.00229687
0.00230417
0.00231146
0.00231875
0.00232604
0.00233333
0.00234062
0.00234792
0.00235521
0.0023625
0.00236979
0.00237708
0.00238437
0.00239167
0.00239896
0.00240625
0.00241354
0.00242083
0.00242812
0.00243542
0.00244271
0.00245
0.00245729
0.00246458
0.00247187
0.00247917
0.00248646
0.00249375
0.00250104
0.00250833
0.00251562
0.00252292
0.00253021
0.0025375
0.00254479
0.00255208
0.00255937
0.00256667
0.00257396
0.00258125
0.00258854
0.00259583
0.00260312
0.00261042
0.00261771
0.002625
0.00263229
0.00263958
0.00264687
0.00265417
0.00266146
0.00266875
0.00267604
0.00268333
0.00269062
0.00269792
0.00270521
0.0027125
0.00271979
0.00272708
0.00273437
0.00274167
0.00274896
0.00275625
0.00276354
0.00277083
0.00277812
0.00278542
0.00279271
0.0028
0.00280729
0.00281458
0.00282187
0.00282917
0.00283646
0.00284375
0.00285104
0.00285833
0.00286562
0.00287292
0.00288021
0.0028875
0.00289479
0.00290208
0.00290937
0.00291667
];

Moment=[0
1.35268e+008
2.07581e+008
2.11132e+008
2.12375e+008
2.1295e+008
2.13263e+008
2.13451e+008
2.13574e+008
2.13658e+008
2.13717e+008
2.13762e+008
2.13796e+008
2.13822e+008
2.13843e+008
2.13859e+008
2.13873e+008
2.13885e+008
2.13894e+008
2.13902e+008
2.13909e+008
2.13915e+008
2.1392e+008
2.13925e+008
2.13929e+008
2.13932e+008
2.13935e+008
2.13938e+008
2.13941e+008
2.1395e+008
2.14027e+008
2.14175e+008
2.14395e+008
2.14695e+008
2.15066e+008
2.15472e+008
2.15881e+008
2.16293e+008
2.16707e+008
2.17124e+008
2.17542e+008
2.17962e+008
2.18383e+008
2.18806e+008
2.1923e+008
2.19655e+008
2.20081e+008
2.20509e+008
2.20937e+008
2.21366e+008
2.21796e+008
2.22226e+008
2.22657e+008
2.23089e+008
2.23521e+008
2.23954e+008
2.24388e+008
2.24821e+008
2.25255e+008
2.2569e+008
2.26125e+008
2.2656e+008
2.26996e+008
2.27431e+008
2.27868e+008
2.28304e+008
2.28741e+008
2.29177e+008
2.29615e+008
2.30052e+008
2.30489e+008
2.30927e+008
2.31365e+008
2.31803e+008
2.32241e+008
2.32679e+008
2.33118e+008
2.33556e+008
2.33995e+008
2.34434e+008
2.34873e+008
2.35312e+008
2.35751e+008
2.3619e+008
2.36629e+008
2.37069e+008
2.37508e+008
2.37948e+008
2.38387e+008
2.38827e+008
2.39267e+008
2.39707e+008
2.40147e+008
2.40587e+008
2.41027e+008
2.41467e+008
2.41907e+008
2.42347e+008
2.42787e+008
2.43228e+008
2.43668e+008
2.44108e+008
2.44549e+008
2.44989e+008
2.4543e+008
2.4587e+008
2.46311e+008
2.46752e+008
2.47192e+008
2.47633e+008
2.48074e+008
2.48514e+008
2.48955e+008
2.49396e+008
2.49837e+008
2.50276e+008
2.50712e+008
2.51143e+008
2.51569e+008
2.51992e+008
2.5241e+008
2.52825e+008
2.53236e+008
2.53643e+008
2.54047e+008
2.54448e+008
2.54844e+008
2.55235e+008
2.55622e+008
2.56005e+008
2.56383e+008
2.56757e+008
2.57127e+008
2.57494e+008
2.57856e+008
2.58215e+008
2.58571e+008
2.58924e+008
2.59278e+008
2.59631e+008
2.59984e+008
2.60336e+008
2.60689e+008
2.61042e+008
2.61394e+008
2.61746e+008
2.62098e+008
2.6245e+008
2.62802e+008
2.63154e+008
2.63506e+008
2.63857e+008
2.64209e+008
2.6456e+008
2.64911e+008
2.65262e+008
2.65613e+008
2.65964e+008
2.66315e+008
2.66665e+008
2.67016e+008
2.67366e+008
2.67717e+008
2.68067e+008
2.68417e+008
2.68768e+008
2.69118e+008
2.69468e+008
2.69818e+008
2.70168e+008
2.70517e+008
2.70867e+008
2.71216e+008
2.71563e+008
2.71908e+008
2.7225e+008
2.72589e+008
2.72927e+008
2.73262e+008
2.73595e+008
2.73925e+008
2.74254e+008
2.7458e+008
2.74905e+008
2.75227e+008
2.75548e+008
2.75867e+008
2.76184e+008
2.76498e+008
2.7681e+008
2.7712e+008
2.77427e+008
2.77731e+008
2.78033e+008
2.78333e+008
2.78631e+008
2.78926e+008
2.79219e+008
2.7951e+008
2.79799e+008
2.80086e+008
2.80371e+008
2.80654e+008
2.80935e+008
2.81214e+008
2.81492e+008
2.81769e+008
2.82047e+008
2.82324e+008
2.82602e+008
2.82879e+008
2.83156e+008
2.83433e+008
2.8371e+008
2.83987e+008
2.84263e+008
2.8454e+008
2.84816e+008
2.85093e+008
2.85369e+008
2.85645e+008
2.85921e+008
2.86197e+008
2.86473e+008
2.86749e+008
2.87024e+008
2.873e+008
2.87575e+008
2.87851e+008
2.88126e+008
2.884e+008
2.88671e+008
2.88941e+008
2.89209e+008
2.89475e+008
2.89739e+008
2.90002e+008
2.90262e+008
2.90521e+008
2.90779e+008
2.91035e+008
2.91289e+008
2.91541e+008
2.91792e+008
2.92042e+008
2.9229e+008
2.92536e+008
2.92781e+008
2.93025e+008
2.93267e+008
2.93508e+008
2.93746e+008
2.93983e+008
2.94218e+008
2.94451e+008
2.94682e+008
2.94911e+008
2.95139e+008
2.95364e+008
2.95588e+008
2.9581e+008
2.96031e+008
2.9625e+008
2.96467e+008
2.96682e+008
2.96896e+008
2.97109e+008
2.97319e+008
2.97529e+008
2.97736e+008
2.97943e+008
2.98148e+008
2.98351e+008
2.98554e+008
2.98756e+008
2.98958e+008
2.9916e+008
2.99363e+008
2.99565e+008
2.99766e+008
2.99968e+008
3.0017e+008
3.00371e+008
3.00573e+008
3.00774e+008
3.00975e+008
3.01176e+008
3.01376e+008
3.01575e+008
3.01771e+008
3.01967e+008
3.0216e+008
3.02353e+008
3.02543e+008
3.02733e+008
3.02921e+008
3.03107e+008
3.03292e+008
3.03476e+008
3.03659e+008
3.0384e+008
3.04019e+008
3.04198e+008
3.04375e+008
3.04551e+008
3.04726e+008
3.04899e+008
3.05071e+008
3.05242e+008
3.05412e+008
3.05581e+008
3.05748e+008
3.05915e+008
3.0608e+008
3.06243e+008
3.06404e+008
3.06565e+008
3.06723e+008
3.0688e+008
3.07036e+008
3.0719e+008
3.07343e+008
3.07494e+008
3.07644e+008
3.07792e+008
3.07939e+008
3.08085e+008
3.08229e+008
3.08372e+008
3.08513e+008
3.08653e+008
3.08792e+008
3.0893e+008
3.09067e+008
3.09202e+008
3.09336e+008
3.09468e+008
3.096e+008
3.0973e+008
3.0986e+008
3.09988e+008
3.10115e+008
3.10242e+008
3.10369e+008
3.10495e+008
3.1062e+008
3.10743e+008
3.10866e+008
3.10986e+008
3.11106e+008
3.11224e+008
3.11342e+008
3.11457e+008
3.11572e+008
3.11685e+008
3.11798e+008
3.11909e+008
3.12019e+008
3.12127e+008
3.12235e+008
3.12341e+008
3.12447e+008
3.12551e+008
3.12654e+008
3.12756e+008
3.12857e+008
3.12957e+008
3.13056e+008
3.13154e+008
3.13251e+008
3.13347e+008
3.13442e+008
3.13536e+008
3.13628e+008
3.1372e+008
3.13811e+008
3.13901e+008
3.13989e+008
3.14075e+008
3.14161e+008
3.14245e+008
3.14328e+008
3.14409e+008
3.1449e+008
3.14569e+008
3.14647e+008
3.14723e+008
3.14799e+008
3.14873e+008
3.14946e+008
3.15018e+008
3.15088e+008
3.15158e+008
3.15226e+008
3.15294e+008
3.1536e+008
3.15425e+008
3.15489e+008
3.15552e+008
3.15614e+008
3.15675e+008
];
%% imaging
figure(1)
IMAGE1=imread('FiberSteelSectionMomentCurvature-image.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE1=imread('FiberIIPlateSteelSectionMomentCurvature.jpg');
image(IMAGE1);axis image;axis off;
%% Plot
figure(3)
P1=plot(ESS,TUSS,[0;ey;esu],[0;fy;fy],'r--');set(P1,'LineWidth',3);
xlabel('Strain');ylabel('Stress (kN/mm^2)')
legend('Exact','Elastic-Perfect Plastic (Simple)','Location','NorthEastOutside');grid on;
title('Steel Section Stress-Strain Diagram','FontSize',12,'color','b')
figure(4)
p1=plot(I,DX,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
figure(5)
p1=plot(I,IT,'b--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
figure(6)
plot(ES,XX,'--','LineWidth',3);xlabel('Top steel strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top steel strain and Neuteral axis ','Color','b');grid on
figure(7)
plot(ES,CUR,'--','LineWidth',2);xlabel('Top steel strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top steel strain and Curvature ','Color','b');grid on
figure(8)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Steel strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the steel fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(9)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Steel stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the steel fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(10)
plot(CUR,EI,'black--','LineWidth',3)
title('#  Flextural Rigidity(EI)-CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(11)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral Axis-Flextural Rigidity(EI) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(12)
plot(Cur,Mom,X,Y,'r--',phi1,Mi1,'g-.','LineWidth',3)
title(['# STEEL SECTION MOMENT CURVATURE DIAGRAM #','  EI_e : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Exact','Exact - bilinear fitted','Elastic-Perfect Plastic (Simple)','Location','NorthEastOutside');grid on;
figure(13)
plot(Cur,Mom,Curvature*1000,Moment*10^-6,'r--','LineWidth',3)
title(['# STEEL SECTION MOMENT CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('MATLAB','C++','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberIIPlateSteelSectionMomentCurvature-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*   Moment-Curvature Analysis of Double I steel sections    *\n');
fprintf(fid,'*   with Double Plates on Flange                            *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*     This program is written by salar delavar ghashghaei   *\n');  
fprintf(fid,'*          E-mail:salar.d.ghashghaei@gmail.com              *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*Unit: Newton-Milimeter                                     *\n');
fprintf(fid,'*Given:Section Properties , Steel Section properties ,      *\n');
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
fprintf(fid,'=Steel Section Moment-Curvature=\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur;Mom]);
fprintf(fid,'+==========================================================+\n');
fprintf(fid,'= Steel Section curve fitted  = Steel Section Simple Anal. =\n');
fprintf(fid,'     Curvature    Moment          Curvature    Moment       \n');
fprintf(fid,'       (1/m)      (kN.m)            (1/m)      (kN.m)       \n');
fprintf(fid,'------------------------------------------------------------\n');
fprintf(fid,'%10.4f %10.3f %20.4f %10.3f\n',[X;Y;phi1';Mi1']);
fprintf(fid,'+==========================================================+\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================================================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I;ES;XX;CUR;EI]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);