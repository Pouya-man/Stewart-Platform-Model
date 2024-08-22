function [ MXp ] = Dynamics( F, thetap,thetab,rp,rb,Px, PPy, Pz, ro, betta, fii, Pxdot, Pydot, Pzdot, rodot, bettadot, fiidot,mpl, Ix, Iy, Iz,g, del, m1, m2,lj1, lj2 )
%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%
clc
%clear all
syms thetap thetab rp rb ro betta fii Px PPy Pz rodot bettadot fiidot Pxdot Pydot Pzdot aldot bedot gadot mpl Ix Iy Iz g del m1 m2 lj1 lj2 k GAxx GAyy GAzz wx wy wz ww Pdot tjj hi vv uix uiy uiz kii F f1 f2 f3 f4 f5 f6 Xpddot alddot beddot gaddot Xddot Yddot Zddot
P = [Px;PPy;Pz];
w = [cos(betta),0,0;0,1,-sin(ro);-sin(betta),0,cos(ro)]*[rodot;bettadot;fiidot];
Pdot=[Pxdot;Pydot;Pzdot];
Xp = [Px;PPy;Pz;ro;betta;fii];
degdot = [rodot;bettadot;fiidot];
Xpdot = [Pxdot;Pydot;Pzdot;rodot;bettadot;fiidot];
%%%%%%%%%%%%%%%%% Inverse kinematics %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% The connection points to the fixed platform (GA(i)) %%%%%%%%%%%%%%%%%
for i=[1,3,5]
landa(i)=((i*(pi/3))+(thetap/2)-(pi/3));
%disp('landa:'), disp(i),disp(landa(i));
end
for i=[2,4,6]
landa(i)=landa(i-1)-thetap+(((2*pi)/3));
%disp('landa:'), disp(i),disp(landa(i))
end
for i=1:6
GAx(i)=rp*cos(landa(i));
GAy(i)=rp*sin(landa(i));
GAz(i)=0;
end
for i=1:6
GA = [GAx(i);GAy(i);GAz(i)];
%disp('GA:'),disp(i),disp(GA)
end
%%%%%%%%%%%%%%%%% The connection points to the moving platform (B(i)) %%%%%%%%%%%%%%%%%
for i=[1,3,5]
delta(i)=((i*(pi/3))-(thetab/2));
%disp('delta:'), disp(i),disp(delta(i))
end
for i=[2,4,6]
delta(i) = delta(i-1)+ thetab;
%disp('delta:'), disp(i),disp(delta(i))
end
for i=1:6
Bx(i)=rb*cos(delta(i));
By(i)=rb*sin(delta(i));
Bz(i)=0;
end
for i=1:6
B=[Bx(i);By(i);Bz(i)];
%disp('B:'),disp(i),disp(B)
end
%%%%%%%%%%%%%%%%%  BRA: rotational matrix %%%%%%%%%%%%%%%%%
Rz=[cos(fii),-sin(fii),0;sin(fii),cos(fii),0;0,0,1];
Ry=[cos(betta),0,sin(betta);0,1,0;-sin(betta),0,cos(betta)];
Rx=[1,0,0;0,cos(ro),-sin(ro);0,sin(ro),cos(ro)];
BRA=Rz*Ry*Rx;
r11=BRA(1,1);
r12=BRA(1,2);
r13=BRA(1,3);
r21=BRA(2,1);
r22=BRA(2,2);
r23=BRA(2,3);
r31=BRA(3,1);
r32=BRA(3,2);
r33=BRA(3,3);
%disp('BRA:'),disp(BRA)
%%%%%%%%%%%%%%%%% END BRA %%%%%%%%%%%%%%%%%
ww=[0,-((cos(ro)*fiidot)+(-sin(betta)*rodot)),(bettadot-(fiidot*sin(ro)));((cos(ro)*fiidot)+(-sin(betta)*rodot)),0,-(cos(betta)*rodot);-(bettadot-(fiidot*sin(ro))),(cos(betta)*rodot),0];
T=(ww*([GAxx;GAyy;GAzz]));
Tjj=Pdot+T;
vv=Tjj'*hi*Tjj;
uii=[uix;uiy;uiz];
Vii=(Pdot'*uii)+(T'*uii);
cc=Vii*kii*Vii;
ll=.5*(m1+m2)*(vv-cc);
%%%%%%%%%%%%%%%%% Li: link vector %%%%%%%%%%%%%%%%%
for i=1:6
L=(BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)];
%disp('L:'),disp(i),disp(L)
end
%%%%%%%%%%%%%%%%% END Li %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  l: length of the arms %%%%%%%%%%%%%%%%%
for i=1:6
l=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
%disp('l:'),disp(i),disp(l)
end
for i=1
l1=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
end
for i=2
l2=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
end
for i=3
l3=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
end
for i=4
l4=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
end
for i=5
l5=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
end
for i=6
l6=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
end
%%%%%%%%%%%%%%%%% end invers kenematic %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START JACOBIN %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% f: unit vector %%%%%%%%%%%%%%%%%
for i=1:6
f=((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5);
%disp('f:'),disp(i),disp(f)
end
for i=1
f1=((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5);
end
for i=2
f2=((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5);
end
for i=3
f3=((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5);
end
for i=4
f4=((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5);
end
for i=5
f5=((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5);
end
for i=6
f6=((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5);
end
%%%%%%%%%%%%%%%%% END u %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START JAI %%%%%%%%%%%%%%%%%
JIA=[f1(1,1),f1(2,1),f1(3,1),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,f2(1,1),f2(2,1),f2(3,1),0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,f3(1,1),f3(2,1),f3(3,1),0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,f4(1,1),f4(2,1),f4(3,1),0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,f5(1,1),f5(2,1),f5(3,1),0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,f6(1,1),f6(2,1),f6(3,1)];
%disp('JIA:'),disp(JIA)
%%%%%%%%%%%%%%%%% END JIA %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START a %%%%%%%%%%%%%%%%%
a=P*(1/(((Px^2)+(PPy^2)+(Pz^2))^.5));
ax=a(1,1);
ay=a(2,1);
az=a(3,1);
S = [0,-az,ay;az,0,-ax;-ay,ax,0];
%disp('S:'),disp(S)
Sx=[0,0,0;0,0,-ax;0,ax,0];
Sy=[0,0,ay;0,0,0;-ay,0,0];
Sz=[0,-az,0;az,0,0;0,0,0];
%%%%%%%%%%%%%%%%% END a and S %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START JIIA %%%%%%%%%%%%%%%%%
GA1=[GAx(1);GAy(1);GAz(1)];
GA2=[GAx(2);GAy(2);GAz(2)];
GA3=[GAx(3);GAy(3);GAz(3)];
GA4=[GAx(4);GAy(4);GAz(4)];
GA5=[GAx(5);GAy(5);GAz(5)];
GA6=[GAx(6);GAy(6);GAz(6)];
s11=((((Ry*Sx)*Rx)*Rz)*GA1);
s12=((((Sy*Ry)*Rx)*Rz)*GA1);
s13=((((Ry*Rx)*Sz)*Rz)*GA1);
s21=((((Ry*Sx)*Rx)*Rz)*GA2);
s22=((((Sy*Ry)*Rx)*Rz)*GA2);
s23=((((Ry*Rx)*Sz)*Rz)*GA2);
s31=((((Ry*Sx)*Rx)*Rz)*GA3);
s32=((((Sy*Ry)*Rx)*Rz)*GA3);
s33=((((Ry*Rx)*Sz)*Rz)*GA3);
s41=((((Ry*Sx)*Rx)*Rz)*GA4);
s42=((((Sy*Ry)*Rx)*Rz)*GA4);
s43=((((Ry*Rx)*Sz)*Rz)*GA4);
s51=((((Ry*Sx)*Rx)*Rz)*GA5);
s52=((((Sy*Ry)*Rx)*Rz)*GA5);
s53=((((Ry*Rx)*Sz)*Rz)*GA5);
s61=((((Ry*Sx)*Rx)*Rz)*GA6);
s62=((((Sy*Ry)*Rx)*Rz)*GA6);
s63=((((Ry*Rx)*Sz)*Rz)*GA6);
JIIA=[1,0,0,s11(1,1),s12(1,1),s13(1,1);0,1,0,s11(2,1),s12(2,1),s13(2,1);0,0,1,s11(3,1),s12(3,1),s13(3,1);1,0,0,s21(1,1),s22(1,1),s23(1,1);0,1,0,s21(2,1),s22(2,1),s23(2,1);0,0,1,s21(3,1),s22(3,1),s23(3,1);1,0,0,s31(1,1),s32(1,1),s33(1,1);0,1,0,s31(2,1),s32(2,1),s33(2,1);0,0,1,s31(3,1),s32(3,1),s33(3,1);1,0,0,s41(1,1),s42(1,1),s43(1,1);0,1,0,s41(2,1),s42(2,1),s43(2,1);0,0,1,s41(3,1),s42(3,1),s43(3,1);1,0,0,s51(1,1),s52(1,1),s53(1,1);0,1,0,s51(2,1),s52(2,1),s53(2,1);0,0,1,s51(3,1),s52(3,1),s53(3,1);1,0,0,s61(1,1),s62(1,1),s63(1,1);0,1,0,s61(2,1),s62(2,1),s63(2,1);0,0,1,s61(3,1),s62(3,1),s63(3,1)];
%disp('JIIA:'),disp(JIIA)
%%%%%%%%%%%%%%%%% END JIIA %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START JA & Lidot1 %%%%%%%%%%%%%%%%%
JA=JIA*JIIA;
Lidot1=JA*Xpdot;
%%%%%%%%%%%%%%%%% END JA & Lidot1 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START JIB %%%%%%%%%%%%%%%%%
k1=cross((BRA*GA1),f1);
k2=cross((BRA*GA2),f2);
k3=cross((BRA*GA3),f3);
k4=cross((BRA*GA4),f4);
k5=cross((BRA*GA5),f5);
k6=cross((BRA*GA6),f6);
JIB =[f1(1,1),f1(2,1),f1(3,1),k1(1,1),k1(2,1),k1(3,1);f2(1,1),f2(2,1),f2(3,1),k2(1,1),k2(2,1),k2(3,1);f3(1,1),f3(2,1),f3(3,1),k3(1,1),k3(2,1),k3(3,1);f4(1,1),f4(2,1),f4(3,1),k4(1,1),k4(2,1),k4(3,1);f5(1,1),f5(2,1),f5(3,1),k5(1,1),k5(2,1),k5(3,1);f6(1,1),f6(2,1),f6(3,1),k6(1,1),k6(2,1),k6(3,1)];
%%%%%%%%%%%%%%%%% END JIB %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START JIIB %%%%%%%%%%%%%%%%%
JIIB =[1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,cos(betta),0,0;0,0,0,0,1,-sin(ro);0,0,0,-sin(betta),0,cos(ro)];
%%%%%%%%%%%%%%%%% END JIIB %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START JA & Lidot1 %%%%%%%%%%%%%%%%%
JB=JIB*JIIB;
Lidot2 = JB*Xpdot;
%%%%%%%%%%%%%%%%% END JA & Lidot1 %%%%%%%%%%%%%%%%%
Lidot1=JA*Xpdot;
%%%%%%%%%%%%%%%%% START Lidot3 %%%%%%%%%%%%%%%%%
for i=1:6
Lidot3 = (dot(Pdot,(((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))))+(dot((cross(w,(BRA*([GAx(i);GAy(i);GAz(i)])))),((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)));
%disp('Lidot3:'),disp(i),disp(Lidot3)
end
for i=1:6
VAj=Pdot+(cross(w,(BRA*[GAx(i);GAy(i);GAz(i)])));
%disp('VAj:'),disp(i),disp(VAj);
end
%%%%%%%%%%%%%%%%% END Lidot3 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% END JACOBIN %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START Dynamics %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Kinetic and potential energies: moving platform %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START Kpl trans %%%%%%%%%%%%%%%%%
Kpltrans = 0.5*(mpl*(Pxdot^2+Pydot^2+Pzdot^2))
%%%%%%%%%%%%%%%%% END Kpl trans %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START Kpl rot %%%%%%%%%%%%%%%%%
Ohmplf = [cos(fii),cos(ro)*sin(fii),(-cos(ro)*cos(fii)*sin(betta))-(cos(ro)*sin(ro)*sin(fii))+(cos(ro)*cos(betta)*sin(ro)*sin(fii)); -sin(fii),cos(ro)*cos(fii),(-cos(ro)*cos(fii)*sin(ro))-(cos(ro)*sin(betta)*sin(fii))+(cos(ro)*cos(betta)*sin(ro)*cos(fii)); 0,-sin(ro),(sin(ro))^2+((cos(ro))^2*cos(betta))]*[rodot;bettadot;fiidot]
Imf = [Ix,0,0;0,Iy,0;0,0,Iz];
Kplrot = 0.5*((Ohmplf'*Imf)*Ohmplf)
Kpl = Kpltrans + Kplrot
disp('Kpl:'),disp(Kpl)
R = [cos(fii),cos(ro)*sin(fii),(-cos(ro)*cos(fii)*sin(betta))-(cos(ro)*sin(ro)*sin(fii))+(cos(ro)*cos(betta)*sin(ro)*sin(fii));-sin(fii),cos(ro)*cos(fii),(-cos(ro)*cos(fii)*sin(ro))-(cos(ro)*sin(betta)*sin(fii))+(cos(ro)*cos(betta)*sin(ro)*cos(fii)); 0,-sin(ro),(sin(ro))^2+((cos(ro))^2*cos(betta))]
RT = R'*Imf*R
Mpl = [mpl 0 0 0 0 0; 0 mpl 0 0 0 0; 0 0 mpl 0 0 0; 0 0 0 RT(1,1) RT(1,2) RT(1,3); 0 0 0 RT(2,1) RT(2,2) RT(2,3); 0 0 0 RT(3,1) RT(3,2) RT(3,3)]
disp('Mpl'),disp(Mpl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Kpl rot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START Ppl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ppl=[0 0 mpl*g 0 0 0]*Xp;
disp('Ppl:'),disp(Ppl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Ppl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% END Kinetic and potential energies: moving platform %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Kinetic energy: legs %%%%%%%%%%%%%%%%%%%%%%%%
Ia = (1/(m1+m2))*((del*m1*lj1)-(0.5*m2*lj2))
for i=1:6
hi=((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2);
disp('hi:'),disp(i),disp(hi)
end
for i=1
h1=((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2);
disp('hi:'),disp(i),disp(hi)
end
for i=2
h2=((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2);
disp('hi:'),disp(i),disp(hi)
end
for i=3
h3=((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2);
disp('hi:'),disp(i),disp(hi)
end
for i=4
h4=((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2);
disp('hi:'),disp(i),disp(hi)
end
for i=5
h5=((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2);
disp('hi:'),disp(i),disp(hi)
end
for i=6
h6=((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2);
disp('hi:'),disp(i),disp(hi)
end
for i=1:6
ki = (((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2
end
for i=1
k1=(((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2
end
for i=2
k2=(((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2
end
for i=3
k3=(((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2
end
for i=4
k4=(((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2
end
for i=5
k5=(((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2
end
for i=6
k6=(((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2
end
Klegs=0;
for i=1:6
Kle= 0.5*(m1+m2)* (((Pdot+(cross(w,(BRA*[GAx(i);GAy(i);GAz(i)]))))' * (((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5)+(m2/(m1+m2)))^2)) * (Pdot+(cross(w,(BRA*[GAx(i);GAy(i);GAz(i)]))))) -(((dot(Pdot,(((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))))+(dot((cross(w,(BRA*([GAx(i);GAy(i);GAz(i)])))),((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))))' * ((((Ia/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))+(m2/(m1+m2)))^2)-(m2/(m1+m2))^2) * ((dot(Pdot,(((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))))+(dot((cross(w,(BRA*([GAx(i);GAy(i);GAz(i)])))),((BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)])/((((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5))))));
Klegs= Klegs + Kle;
disp('Kle:'),disp(i),disp(Kle)
end
MPxdot=(h1-(k1*(f1(1,1))^2))+(h2-(k2*(f2(1,1))^2))+(h3-(k3*(f3(1,1))^2))+(h4-(k4*(f4(1,1))^2))+(h5-(k5*(f5(1,1))^2))+(h6-(k6*(f6(1,1))^2));
MPydot=(h1-(k1*(f1(2,1))^2))+(h2-(k2*(f2(2,1))^2))+(h3-(k3*(f3(2,1))^2))+(h4-(k4*(f4(2,1))^2))+(h5-(k5*(f5(2,1))^2))+(h6-(k6*(f6(2,1))^2));
MPzdot=(h1-(k1*(f1(3,1))^2))+(h2-(k2*(f2(3,1))^2))+(h3-(k3*(f3(3,1))^2))+(h4-(k4*(f4(3,1))^2))+(h5-(k5*(f5(3,1))^2))+(h6-(k6*(f6(3,1))^2));
Maldot=((h1*sin(betta)^2*GAy(1)^2)+(h1*sin(betta)^2*GAx(1)^2)-(h1*cos(betta)^2*GAz(1)^2)+(h1*cos(betta)^2*GAy(1)^2)-(k1*sin(betta)^2*f1(1,1)^2*GAy(1)^2)-(k1*sin(betta)^2*f1(2,1)^2*GAx(1)^2)-(k1*cos(betta)^2*f1(2,1)^2*GAz(1)^2)-(k1*cos(betta)^2*f1(3,1)^2*GAy(1)^2)+(2*h1*sin(betta)*GAx(1)*cos(betta)*GAz(1))+(2*k1*sin(betta)^2*f1(1,1)*f1(2,1)*GAy(1)*GAx(1))+(2*k1*sin(betta)*cos(betta)*f1(1,1)*f1(2,1)*GAy(1)*GAz(1))-(2*k1*sin(betta)*cos(betta)*f1(1,1)*f1(3,1)*GAy(1)*GAy(1))-(2*k1*sin(betta)*cos(betta)*f1(2,1)*f1(2,1)*GAx(1)*GAz(1))+(2*k1*sin(betta)*cos(betta)*f1(2,1)*f1(3,1)*GAx(1)*GAy(1))+(2*k1*cos(betta)*cos(betta)*f1(2,1)*f1(3,1)*GAz(1)*GAy(1))+(2*k1*cos(betta)*cos(betta)*f1(2,1)*f1(3,1)*GAz(1)*GAy(1)))+((h2*sin(betta)^2*GAy(2)^2)+(h2*sin(betta)^2*GAx(2)^2)-(h2*cos(betta)^2*GAz(2)^2)+(h2*cos(betta)^2*GAy(2)^2)-(k2*sin(betta)^2*f2(1,1)^2*GAy(2)^2)-(k2*sin(betta)^2*f2(2,1)^2*GAx(2)^2)-(k2*cos(betta)^2*f2(2,1)^2*GAz(2)^2)-(k2*cos(betta)^2*f2(3,1)^2*GAy(2)^2)+(2*h2*sin(betta)*GAx(2)*cos(betta)*GAz(2))+(2*k2*sin(betta)^2*f2(1,1)*f2(2,1)*GAy(2)*GAx(2))+(2*k2*sin(betta)*cos(betta)*f2(1,1)*f2(2,1)*GAy(2)*GAz(2))-(2*k2*sin(betta)*cos(betta)*f2(1,1)*f2(3,1)*GAy(2)*GAy(2))-(2*k2*sin(betta)*cos(betta)*f2(2,1)*f2(2,1)*GAx(2)*GAz(2))+(2*k2*sin(betta)*cos(betta)*f2(2,1)*f2(3,1)*GAx(2)*GAy(2))+(2*k2*cos(betta)*cos(betta)*f2(2,1)*f2(3,1)*GAz(2)*GAy(2))+(2*k2*cos(betta)*cos(betta)*f2(2,1)*f2(3,1)*GAz(2)*GAy(2)))+((h3*sin(betta)^2*GAy(3)^2)+(h3*sin(betta)^2*GAx(3)^2)-(h3*cos(betta)^2*GAz(3)^2)+(h3*cos(betta)^2*GAy(3)^2)-(k3*sin(betta)^2*f3(1,1)^2*GAy(3)^2)-(k3*sin(betta)^2*f3(2,1)^2*GAx(3)^2)-(k3*cos(betta)^2*f3(2,1)^2*GAz(3)^2)-(k3*cos(betta)^2*f3(3,1)^2*GAy(3)^2)+(2*h3*sin(betta)*GAx(3)*cos(betta)*GAz(3))+(2*k3*sin(betta)^2*f3(1,1)*f3(2,1)*GAy(3)*GAx(3))+(2*k3*sin(betta)*cos(betta)*f3(1,1)*f3(2,1)*GAy(3)*GAz(3))-(2*k3*sin(betta)*cos(betta)*f3(1,1)*f3(3,1)*GAy(3)*GAy(3))-(2*k3*sin(betta)*cos(betta)*f3(2,1)*f3(2,1)*GAx(3)*GAz(3))+(2*k3*sin(betta)*cos(betta)*f3(2,1)*f3(3,1)*GAx(3)*GAy(3))+(2*k3*cos(betta)*cos(betta)*f3(2,1)*f3(3,1)*GAz(3)*GAy(3))+(2*k3*cos(betta)*cos(betta)*f3(2,1)*f3(3,1)*GAz(3)*GAy(3)))+((h4*sin(betta)^2*GAy(4)^2)+(h4*sin(betta)^2*GAx(4)^2)-(h4*cos(betta)^2*GAz(4)^2)+(h4*cos(betta)^2*GAy(4)^2)-(k4*sin(betta)^2*f4(1,1)^2*GAy(4)^2)-(k4*sin(betta)^2*f4(2,1)^2*GAx(4)^2)-(k4*cos(betta)^2*f4(2,1)^2*GAz(4)^2)-(k4*cos(betta)^2*f4(3,1)^2*GAy(4)^2)+(2*h4*sin(betta)*GAx(4)*cos(betta)*GAz(4))+(2*k4*sin(betta)^2*f4(1,1)*f4(2,1)*GAy(4)*GAx(4))+(2*k4*sin(betta)*cos(betta)*f4(1,1)*f4(2,1)*GAy(4)*GAz(4))-(2*k4*sin(betta)*cos(betta)*f4(1,1)*f4(3,1)*GAy(4)*GAy(4))-(2*k4*sin(betta)*cos(betta)*f4(2,1)*f4(2,1)*GAx(4)*GAz(4))+(2*k4*sin(betta)*cos(betta)*f4(2,1)*f4(3,1)*GAx(4)*GAy(4))+(2*k4*cos(betta)*cos(betta)*f4(2,1)*f4(3,1)*GAz(4)*GAy(4))+(2*k4*cos(betta)*cos(betta)*f4(2,1)*f4(3,1)*GAz(4)*GAy(4)))+((h5*sin(betta)^2*GAy(5)^2)+(h5*sin(betta)^2*GAx(5)^2)-(h5*cos(betta)^2*GAz(5)^2)+(h5*cos(betta)^2*GAy(5)^2)-(k5*sin(betta)^2*f5(1,1)^2*GAy(5)^2)-(k5*sin(betta)^2*f5(2,1)^2*GAx(5)^2)-(k5*cos(betta)^2*f5(2,1)^2*GAz(5)^2)-(k5*cos(betta)^2*f5(3,1)^2*GAy(5)^2)+(2*h5*sin(betta)*GAx(5)*cos(betta)*GAz(5))+(2*k5*sin(betta)^2*f5(1,1)*f5(2,1)*GAy(5)*GAx(5))+(2*k5*sin(betta)*cos(betta)*f5(1,1)*f5(2,1)*GAy(5)*GAz(5))-(2*k5*sin(betta)*cos(betta)*f5(1,1)*f5(3,1)*GAy(5)*GAy(5))-(2*k5*sin(betta)*cos(betta)*f5(2,1)*f5(2,1)*GAx(5)*GAz(5))+(2*k5*sin(betta)*cos(betta)*f5(2,1)*f5(3,1)*GAx(5)*GAy(5))+(2*k5*cos(betta)*cos(betta)*f5(2,1)*f5(3,1)*GAz(5)*GAy(5))+(2*k5*cos(betta)*cos(betta)*f5(2,1)*f5(3,1)*GAz(5)*GAy(5)))+((h6*sin(betta)^2*GAy(6)^2)+(h6*sin(betta)^2*GAx(6)^2)-(h6*cos(betta)^2*GAz(6)^2)+(h6*cos(betta)^2*GAy(6)^2)-(k6*sin(betta)^2*f6(1,1)^2*GAy(6)^2)-(k6*sin(betta)^2*f6(2,1)^2*GAx(6)^2)-(k6*cos(betta)^2*f6(2,1)^2*GAz(6)^2)-(k6*cos(betta)^2*f6(3,1)^2*GAy(6)^2)+(2*h6*sin(betta)*GAx(6)*cos(betta)*GAz(6))+(2*k6*sin(betta)^2*f6(1,1)*f6(2,1)*GAy(6)*GAx(6))+(2*k6*sin(betta)*cos(betta)*f6(1,1)*f6(2,1)*GAy(6)*GAz(6))-(2*k6*sin(betta)*cos(betta)*f6(1,1)*f6(3,1)*GAy(6)*GAy(6))-(2*k6*sin(betta)*cos(betta)*f6(2,1)*f6(2,1)*GAx(6)*GAz(6))+(2*k6*sin(betta)*cos(betta)*f6(2,1)*f6(3,1)*GAx(6)*GAy(6))+(2*k6*cos(betta)*cos(betta)*f6(2,1)*f6(3,1)*GAz(6)*GAy(6))+(2*k6*cos(betta)*cos(betta)*f6(2,1)*f6(3,1)*GAz(6)*GAy(6)));
Mbedot=(h1*GAz(1)^2)+(h1*GAx(1)^2)+(h1*GAz(1)^2*f1(1,1)^2)+(h1*GAx(1)^2*f1(3,1)^2)+(2*k1*GAz(1)*GAx(1)*f1(1,1)*f1(3,1))+(h2*GAz(2)^2)+(h2*GAx(2)^2)+(h2*GAz(2)^2*f2(1,1)^2)+(h2*GAx(2)^2*f2(3,1)^2)+(2*k2*GAz(2)*GAx(2)*f2(1,1)*f2(3,1))+(h3*GAz(3)^2)+(h3*GAx(3)^2)+(h3*GAz(3)^2*f3(1,1)^2)+(h3*GAx(3)^2*f3(3,1)^2)+(2*k3*GAz(3)*GAx(3)*f3(1,1)*f3(3,1))+(h4*GAz(4)^2)+(h4*GAx(4)^2)+(h4*GAz(4)^2*f4(1,1)^2)+(h4*GAx(4)^2*f4(3,1)^2)+(2*k4*GAz(4)*GAx(4)*f4(1,1)*f4(3,1))+(h5*GAz(5)^2)+(h5*GAx(5)^2)+(h5*GAz(5)^2*f5(1,1)^2)+(h5*GAx(5)^2*f5(3,1)^2)+(2*k5*GAz(5)*GAx(5)*f5(1,1)*f5(3,1))+(h6*GAz(6)^2)+(h6*GAx(6)^2)+(h6*GAz(6)^2*f6(1,1)^2)+(h6*GAx(1)^2*f6(3,1)^2)+(2*k6*GAz(6)*GAx(6)*f6(1,1)*f6(3,1));
Mgadot=(h1*cos(ro)^2*GAy(1)^2)+(h1*sin(ro)^2*GAz(1)^2)+(h1*cos(ro)^2*GAx(1)^2)+(h1*sin(ro)^2*GAx(1)^2)-(k1*cos(ro)^2*GAy(1)^2*f1(1,1)^2)-(k1*sin(ro)^2*GAx(1)^2*f1(1,1)^2)-(k1*cos(ro)^2*GAx(1)^2*f1(2,1)^2)-(k1*sin(ro)^2*GAx(1)^2*f1(3,1)^2)+(2*h1*cos(ro)*sin(ro)*GAy(1)*GAz(1))-(2*k1*cos(ro)*sin(ro)*f1(1,1)^2*GAy(1)*GAx(1))+(2*k1*cos(ro)*cos(ro)*f1(1,1)*f1(2,1)*GAy(1)*GAx(1))+(2*k1*cos(ro)*sin(ro)*f1(1,1)*f1(3,1)*GAy(1)*GAx(1))+(2*k1*cos(ro)*sin(ro)*f1(1,1)*f1(2,1)*GAx(1)*GAx(1))+(2*k1*sin(ro)*sin(ro)*f1(1,1)*f1(3,1)*GAx(1)*GAx(1))-(2*k1*cos(ro)*sin(ro)*f1(3,1)*f1(2,1)*GAx(1)*GAx(1))+(h2*cos(ro)^2*GAy(2)^2)+(h2*sin(ro)^2*GAz(2)^2)+(h2*cos(ro)^2*GAx(2)^2)+(h2*sin(ro)^2*GAx(2)^2)-(k2*cos(ro)^2*GAy(2)^2*f2(1,1)^2)-(k2*sin(ro)^2*GAx(2)^2*f2(1,1)^2)-(k2*cos(ro)^2*GAx(2)^2*f2(2,1)^2)-(k2*sin(ro)^2*GAx(2)^2*f2(3,1)^2)+(2*h2*cos(ro)*sin(ro)*GAy(2)*GAz(2))-(2*k2*cos(ro)*sin(ro)*f2(1,1)^2*GAy(2)*GAx(2))+(2*k2*cos(ro)*cos(ro)*f2(1,1)*f2(2,1)*GAy(2)*GAx(2))+(2*k2*cos(ro)*sin(ro)*f2(1,1)*f2(3,1)*GAy(2)*GAx(2))+(2*k2*cos(ro)*sin(ro)*f2(1,1)*f2(2,1)*GAx(2)*GAx(2))+(2*k2*sin(ro)*sin(ro)*f2(1,1)*f2(3,1)*GAx(2)*GAx(2))-(2*k2*cos(ro)*sin(ro)*f2(3,1)*f2(2,1)*GAx(2)*GAx(2))+(h3*cos(ro)^2*GAy(3)^2)+(h3*sin(ro)^2*GAz(3)^2)+(h3*cos(ro)^2*GAx(3)^2)+(h3*sin(ro)^2*GAx(3)^2)-(k3*cos(ro)^2*GAy(3)^2*f3(1,1)^2)-(k3*sin(ro)^2*GAx(3)^2*f3(1,1)^2)-(k3*cos(ro)^2*GAx(3)^2*f3(2,1)^2)-(k3*sin(ro)^2*GAx(3)^2*f3(3,1)^2)+(2*h3*cos(ro)*sin(ro)*GAy(3)*GAz(3))-(2*k3*cos(ro)*sin(ro)*f3(1,1)^2*GAy(3)*GAx(3))+(2*k3*cos(ro)*cos(ro)*f3(1,1)*f3(2,1)*GAy(3)*GAx(3))+(2*k3*cos(ro)*sin(ro)*f3(1,1)*f3(3,1)*GAy(3)*GAx(3))+(2*k3*cos(ro)*sin(ro)*f3(1,1)*f3(2,1)*GAx(3)*GAx(3))+(2*k3*sin(ro)*sin(ro)*f3(1,1)*f3(3,1)*GAx(3)*GAx(3))-(2*k3*cos(ro)*sin(ro)*f3(3,1)*f3(2,1)*GAx(3)*GAx(3))+(h4*cos(ro)^2*GAy(4)^2)+(h4*sin(ro)^2*GAz(4)^2)+(h4*cos(ro)^2*GAx(4)^2)+(h4*sin(ro)^2*GAx(4)^2)-(k4*cos(ro)^2*GAy(4)^2*f4(1,1)^2)-(k4*sin(ro)^2*GAx(4)^2*f4(1,1)^2)-(k4*cos(ro)^2*GAx(4)^2*f4(2,1)^2)-(k4*sin(ro)^2*GAx(4)^2*f4(3,1)^2)+(2*h4*cos(ro)*sin(ro)*GAy(4)*GAz(4))-(2*k4*cos(ro)*sin(ro)*f4(1,1)^2*GAy(4)*GAx(4))+(2*k4*cos(ro)*cos(ro)*f4(1,1)*f4(2,1)*GAy(4)*GAx(4))+(2*k4*cos(ro)*sin(ro)*f4(1,1)*f4(3,1)*GAy(4)*GAx(4))+(2*k4*cos(ro)*sin(ro)*f4(1,1)*f4(2,1)*GAx(4)*GAx(4))+(2*k4*sin(ro)*sin(ro)*f4(1,1)*f4(3,1)*GAx(4)*GAx(4))-(2*k4*cos(ro)*sin(ro)*f4(3,1)*f4(2,1)*GAx(4)*GAx(4))+(h5*cos(ro)^2*GAy(5)^2)+(h5*sin(ro)^2*GAz(5)^2)+(h5*cos(ro)^2*GAx(5)^2)+(h5*sin(ro)^2*GAx(5)^2)-(k5*cos(ro)^2*GAy(5)^2*f5(1,1)^2)-(k5*sin(ro)^2*GAx(5)^2*f5(1,1)^2)-(k5*cos(ro)^2*GAx(5)^2*f5(2,1)^2)-(k5*sin(ro)^2*GAx(5)^2*f5(3,1)^2)+(2*h5*cos(ro)*sin(ro)*GAy(5)*GAz(5))-(2*k5*cos(ro)*sin(ro)*f5(1,1)^2*GAy(5)*GAx(5))+(2*k5*cos(ro)*cos(ro)*f5(1,1)*f5(2,1)*GAy(5)*GAx(5))+(2*k5*cos(ro)*sin(ro)*f5(1,1)*f5(3,1)*GAy(5)*GAx(5))+(2*k5*cos(ro)*sin(ro)*f5(1,1)*f5(2,1)*GAx(5)*GAx(5))+(2*k5*sin(ro)*sin(ro)*f5(1,1)*f5(3,1)*GAx(5)*GAx(5))-(2*k5*cos(ro)*sin(ro)*f5(3,1)*f5(2,1)*GAx(5)*GAx(5))+(h6*cos(ro)^2*GAy(6)^2)+(h6*sin(ro)^2*GAz(6)^2)+(h6*cos(ro)^2*GAx(6)^2)+(h6*sin(ro)^2*GAx(6)^2)-(k6*cos(ro)^2*GAy(6)^2*f6(1,1)^2)-(k6*sin(ro)^2*GAx(6)^2*f6(1,1)^2)-(k6*cos(ro)^2*GAx(6)^2*f6(2,1)^2)-(k6*sin(ro)^2*GAx(6)^2*f6(3,1)^2)+(2*h6*cos(ro)*sin(ro)*GAy(6)*GAz(6))-(2*k6*cos(ro)*sin(ro)*f6(1,1)^2*GAy(6)*GAx(6))+(2*k6*cos(ro)*cos(ro)*f6(1,1)*f6(2,1)*GAy(6)*GAx(6))+(2*k6*cos(ro)*sin(ro)*f6(1,1)*f6(3,1)*GAy(6)*GAx(6))+(2*k6*cos(ro)*sin(ro)*f6(1,1)*f6(2,1)*GAx(6)*GAx(6))+(2*k6*sin(ro)*sin(ro)*f6(1,1)*f6(3,1)*GAx(6)*GAx(6))-(2*k6*cos(ro)*sin(ro)*f6(3,1)*f6(2,1)*GAx(6)*GAx(6));
Malbedot=(2*h1*sin(betta)*GAy(1)*GAz(1))-(2*h1*cos(betta)*GAy(1)*GAx(1))-(2*k1*sin(betta)*f1(1,1)*f1(1,1)*GAy(1)*GAz(1))+(2*k1*sin(betta)*f1(1,1)*f1(3,1)*GAy(1)*GAx(1))+(2*k1*sin(betta)*f1(1,1)*f1(2,1)*GAx(1)*GAz(1))-(2*k1*cos(betta)*f1(1,1)*f1(3,1)*GAy(1)*GAz(1))+(2*k1*cos(betta)*f1(1,1)*f1(3,1)*GAz(1)*GAz(1))-(2*k1*sin(betta)*f1(2,1)*f1(3,1)*GAx(1)*GAx(1))-(2*k1*cos(betta)*f1(3,1)*f1(2,1)*GAx(1)*GAz(1))+(2*k1*cos(betta)*f1(3,1)*f1(3,1)*GAx(1)*GAy(1))+(2*h2*sin(betta)*GAy(2)*GAz(2))-(2*h2*cos(betta)*GAy(2)*GAx(2))-(2*k2*sin(betta)*f2(1,1)*f2(1,1)*GAy(2)*GAz(2))+(2*k2*sin(betta)*f2(1,1)*f2(3,1)*GAy(2)*GAx(2))+(2*k2*sin(betta)*f2(1,1)*f2(2,1)*GAx(2)*GAz(2))-(2*k2*cos(betta)*f2(1,1)*f2(3,1)*GAy(2)*GAz(2))+(2*k2*cos(betta)*f2(1,1)*f2(3,1)*GAz(2)*GAz(2))-(2*k2*sin(betta)*f2(2,1)*f2(3,1)*GAx(2)*GAx(2))-(2*k2*cos(betta)*f2(3,1)*f2(2,1)*GAx(2)*GAz(2))+(2*k2*cos(betta)*f2(3,1)*f2(3,1)*GAx(2)*GAy(2))+(2*h3*sin(betta)*GAy(3)*GAz(3))-(2*h3*cos(betta)*GAy(3)*GAx(3))-(2*k3*sin(betta)*f3(1,1)*f3(1,1)*GAy(3)*GAz(3))+(2*k3*sin(betta)*f3(1,1)*f3(3,1)*GAy(3)*GAx(3))+(2*k3*sin(betta)*f3(1,1)*f3(2,1)*GAx(3)*GAz(3))-(2*k3*cos(betta)*f3(1,1)*f3(3,1)*GAy(3)*GAz(3))+(2*k3*cos(betta)*f3(1,1)*f3(3,1)*GAz(3)*GAz(3))-(2*k3*sin(betta)*f3(2,1)*f3(3,1)*GAx(3)*GAx(3))-(2*k3*cos(betta)*f3(3,1)*f3(2,1)*GAx(3)*GAz(3))+(2*k3*cos(betta)*f3(3,1)*f3(3,1)*GAx(3)*GAy(3))+(2*h4*sin(betta)*GAy(4)*GAz(4))-(2*h4*cos(betta)*GAy(4)*GAx(4))-(2*k4*sin(betta)*f4(1,1)*f4(1,1)*GAy(4)*GAz(4))+(2*k4*sin(betta)*f4(1,1)*f4(3,1)*GAy(4)*GAx(4))+(2*k4*sin(betta)*f4(1,1)*f4(2,1)*GAx(4)*GAz(4))-(2*k4*cos(betta)*f4(1,1)*f4(3,1)*GAy(4)*GAz(4))+(2*k4*cos(betta)*f4(1,1)*f4(3,1)*GAz(4)*GAz(4))-(2*k4*sin(betta)*f4(2,1)*f4(3,1)*GAx(4)*GAx(4))-(2*k4*cos(betta)*f4(3,1)*f4(2,1)*GAx(4)*GAz(4))+(2*k4*cos(betta)*f4(3,1)*f4(3,1)*GAx(4)*GAy(4))+(2*h5*sin(betta)*GAy(5)*GAz(5))-(2*h5*cos(betta)*GAy(5)*GAx(5))-(2*k5*sin(betta)*f5(1,1)*f5(1,1)*GAy(5)*GAz(5))+(2*k5*sin(betta)*f5(1,1)*f5(3,1)*GAy(5)*GAx(5))+(2*k5*sin(betta)*f5(1,1)*f5(2,1)*GAx(5)*GAz(5))-(2*k5*cos(betta)*f5(1,1)*f5(3,1)*GAy(5)*GAz(5))+(2*k5*cos(betta)*f5(1,1)*f5(3,1)*GAz(5)*GAz(5))-(2*k5*sin(betta)*f5(2,1)*f5(3,1)*GAx(5)*GAx(5))-(2*k5*cos(betta)*f5(3,1)*f5(2,1)*GAx(5)*GAz(5))+(2*k5*cos(betta)*f5(3,1)*f5(3,1)*GAx(5)*GAy(5))+(2*h6*sin(betta)*GAy(6)*GAz(6))-(2*h6*cos(betta)*GAy(6)*GAx(6))-(2*k6*sin(betta)*f6(1,1)*f6(1,1)*GAy(6)*GAz(6))+(2*k6*sin(betta)*f6(1,1)*f6(3,1)*GAy(6)*GAx(6))+(2*k6*sin(betta)*f6(1,1)*f6(2,1)*GAx(6)*GAz(6))-(2*k6*cos(betta)*f6(1,1)*f6(3,1)*GAy(6)*GAz(6))+(2*k6*cos(betta)*f6(1,1)*f6(3,1)*GAz(6)*GAz(6))-(2*k6*sin(betta)*f6(2,1)*f6(3,1)*GAx(6)*GAx(6))-(2*k6*cos(betta)*f6(3,1)*f6(2,1)*GAx(6)*GAz(6))+(2*k6*cos(betta)*f6(3,1)*f6(3,1)*GAx(6)*GAy(6));
Malgadot=(-2*h1*sin(betta)*cos(ro)*GAy(1)*GAy(1))-(2*h1*sin(betta)*sin(ro)*GAy(1)*GAz(1))-(2*h1*cos(betta)*cos(ro)*GAx(1)*GAz(1))+(2*h1*cos(betta)*sin(ro)*GAy(1)*GAx(1))+(2*k1*sin(betta)*sin(ro)*f1(1,1)*f1(1,1)*GAy(1)*GAx(1))-(2*k1*sin(betta)*sin(ro)*f1(1,1)*f1(3,1)*GAy(1)*GAx(1))-(2*k1*sin(betta)*cos(ro)*f1(1,1)*f1(2,1)*GAy(1)*GAx(1))+(2*k1*cos(betta)*cos(ro)*f1(1,1)*f1(3,1)*GAy(1)*GAy(1))-(2*k1*sin(betta)*cos(ro)*f1(1,1)*f1(2,1)*GAy(1)*GAx(1))-(2*k1*cos(betta)*cos(ro)*f1(1,1)*f1(2,1)*GAy(1)*GAz(1))-(2*k1*sin(betta)*sin(ro)*f1(1,1)*f1(2,1)*GAx(1)*GAx(1))+(2*k1*cos(betta)*sin(ro)*f1(1,1)*f1(3,1)*GAy(1)*GAx(1))+(2*k1*sin(betta)*cos(ro)*f1(2,1)*f1(2,1)*GAx(1)*GAx(1))+(2*k1*sin(betta)*sin(ro)*f1(2,1)*f1(3,1)*GAx(1)*GAx(1))-(2*k1*cos(betta)*cos(ro)*f1(2,1)*f1(3,1)*GAx(1)*GAy(1))+(2*k1*cos(betta)*cos(ro)*f1(2,1)*f1(2,1)*GAz(1)*GAx(1))+(2*k1*cos(betta)*sin(ro)*f1(3,1)*f1(2,1)*GAz(1)*GAx(1))-(2*k1*cos(betta)*sin(ro)*f1(1,1)*f1(2,1)*GAz(1)*GAx(1))-(2*k1*cos(betta)*sin(ro)*f1(3,1)*f1(3,1)*GAy(1)*GAx(1))+(2*k1*sin(betta)*cos(ro)*f1(1,1)*f1(1,1)*GAy(1)*GAy(1))-(2*h1*sin(betta)*cos(ro)*GAx(1)*GAx(1))+(-2*h2*sin(betta)*cos(ro)*GAy(2)*GAy(2))-(2*h2*sin(betta)*sin(ro)*GAy(2)*GAz(2))-(2*h2*cos(betta)*cos(ro)*GAx(2)*GAz(2))+(2*h2*cos(betta)*sin(ro)*GAy(2)*GAx(2))+(2*k2*sin(betta)*sin(ro)*f2(1,1)*f2(1,1)*GAy(2)*GAx(2))-(2*k2*sin(betta)*sin(ro)*f2(1,1)*f2(3,1)*GAy(2)*GAx(2))-(2*k2*sin(betta)*cos(ro)*f2(1,1)*f2(2,1)*GAy(2)*GAx(2))+(2*k2*cos(betta)*cos(ro)*f2(1,1)*f2(3,1)*GAy(2)*GAy(2))-(2*k2*sin(betta)*cos(ro)*f2(1,1)*f2(2,1)*GAy(2)*GAx(2))-(2*k2*cos(betta)*cos(ro)*f2(1,1)*f2(2,1)*GAy(2)*GAz(2))-(2*k2*sin(betta)*sin(ro)*f2(1,1)*f2(2,1)*GAx(2)*GAx(2))+(2*k2*cos(betta)*sin(ro)*f2(1,1)*f2(3,1)*GAy(2)*GAx(2))+(2*k2*sin(betta)*cos(ro)*f2(2,1)*f2(2,1)*GAx(2)*GAx(2))+(2*k2*sin(betta)*sin(ro)*f2(2,1)*f2(3,1)*GAx(2)*GAx(2))-(2*k2*cos(betta)*cos(ro)*f2(2,1)*f2(3,1)*GAx(2)*GAy(2))+(2*k2*cos(betta)*cos(ro)*f2(2,1)*f2(2,1)*GAz(2)*GAx(2))+(2*k2*cos(betta)*sin(ro)*f2(3,1)*f2(2,1)*GAz(2)*GAx(2))-(2*k2*cos(betta)*sin(ro)*f2(1,1)*f2(2,1)*GAz(2)*GAx(2))-(2*k2*cos(betta)*sin(ro)*f2(3,1)*f2(3,1)*GAy(2)*GAx(2))+(2*k2*sin(betta)*cos(ro)*f2(1,1)*f2(1,1)*GAy(2)*GAy(2))-(2*h2*sin(betta)*cos(ro)*GAx(2)*GAx(2))+(-2*h3*sin(betta)*cos(ro)*GAy(3)*GAy(3))-(2*h3*sin(betta)*sin(ro)*GAy(3)*GAz(3))-(2*h3*cos(betta)*cos(ro)*GAx(3)*GAz(3))+(2*h3*cos(betta)*sin(ro)*GAy(3)*GAx(3))+(2*k3*sin(betta)*sin(ro)*f3(1,1)*f3(1,1)*GAy(3)*GAx(3))-(2*k3*sin(betta)*sin(ro)*f3(1,1)*f3(3,1)*GAy(3)*GAx(3))-(2*k3*sin(betta)*cos(ro)*f3(1,1)*f3(2,1)*GAy(3)*GAx(3))+(2*k3*cos(betta)*cos(ro)*f3(1,1)*f3(3,1)*GAy(3)*GAy(3))-(2*k3*sin(betta)*cos(ro)*f3(1,1)*f3(2,1)*GAy(3)*GAx(3))-(2*k3*cos(betta)*cos(ro)*f3(1,1)*f3(2,1)*GAy(3)*GAz(3))-(2*k3*sin(betta)*sin(ro)*f3(1,1)*f3(2,1)*GAx(3)*GAx(3))+(2*k3*cos(betta)*sin(ro)*f3(1,1)*f3(3,1)*GAy(3)*GAx(3))+(2*k3*sin(betta)*cos(ro)*f3(2,1)*f3(2,1)*GAx(3)*GAx(3))+(2*k3*sin(betta)*sin(ro)*f3(2,1)*f3(3,1)*GAx(3)*GAx(3))-(2*k3*cos(betta)*cos(ro)*f3(2,1)*f3(3,1)*GAx(3)*GAy(3))+(2*k3*cos(betta)*cos(ro)*f3(2,1)*f3(2,1)*GAz(3)*GAx(3))+(2*k3*cos(betta)*sin(ro)*f3(3,1)*f3(2,1)*GAz(3)*GAx(3))-(2*k3*cos(betta)*sin(ro)*f3(1,1)*f3(2,1)*GAz(3)*GAx(3))-(2*k3*cos(betta)*sin(ro)*f3(3,1)*f3(3,1)*GAy(3)*GAx(3))+(2*k3*sin(betta)*cos(ro)*f3(1,1)*f3(1,1)*GAy(3)*GAy(3))-(2*h3*sin(betta)*cos(ro)*GAx(3)*GAx(3))+(-2*h4*sin(betta)*cos(ro)*GAy(4)*GAy(4))-(2*h4*sin(betta)*sin(ro)*GAy(4)*GAz(4))-(2*h4*cos(betta)*cos(ro)*GAx(4)*GAz(4))+(2*h4*cos(betta)*sin(ro)*GAy(4)*GAx(4))+(2*k4*sin(betta)*sin(ro)*f4(1,1)*f4(1,1)*GAy(4)*GAx(4))-(2*k4*sin(betta)*sin(ro)*f4(1,1)*f4(3,1)*GAy(4)*GAx(4))-(2*k4*sin(betta)*cos(ro)*f4(1,1)*f4(2,1)*GAy(4)*GAx(4))+(2*k4*cos(betta)*cos(ro)*f4(1,1)*f4(3,1)*GAy(4)*GAy(4))-(2*k4*sin(betta)*cos(ro)*f4(1,1)*f4(2,1)*GAy(4)*GAx(4))-(2*k4*cos(betta)*cos(ro)*f4(1,1)*f4(2,1)*GAy(4)*GAz(4))-(2*k4*sin(betta)*sin(ro)*f4(1,1)*f4(2,1)*GAx(4)*GAx(4))+(2*k4*cos(betta)*sin(ro)*f4(1,1)*f4(3,1)*GAy(4)*GAx(4))+(2*k4*sin(betta)*cos(ro)*f4(2,1)*f4(2,1)*GAx(4)*GAx(4))+(2*k4*sin(betta)*sin(ro)*f4(2,1)*f4(3,1)*GAx(4)*GAx(4))-(2*k4*cos(betta)*cos(ro)*f4(2,1)*f4(3,1)*GAx(4)*GAy(4))+(2*k4*cos(betta)*cos(ro)*f4(2,1)*f4(2,1)*GAz(4)*GAx(4))+(2*k4*cos(betta)*sin(ro)*f4(3,1)*f4(2,1)*GAz(4)*GAx(4))-(2*k4*cos(betta)*sin(ro)*f4(1,1)*f4(2,1)*GAz(4)*GAx(4))-(2*k4*cos(betta)*sin(ro)*f4(3,1)*f4(3,1)*GAy(4)*GAx(4))+(2*k4*sin(betta)*cos(ro)*f4(1,1)*f4(1,1)*GAy(4)*GAy(4))-(2*h4*sin(betta)*cos(ro)*GAx(4)*GAx(4))+(-2*h5*sin(betta)*cos(ro)*GAy(5)*GAy(5))-(2*h5*sin(betta)*sin(ro)*GAy(5)*GAz(5))-(2*h5*cos(betta)*cos(ro)*GAx(5)*GAz(5))+(2*h5*cos(betta)*sin(ro)*GAy(5)*GAx(5))+(2*k5*sin(betta)*sin(ro)*f5(1,1)*f5(1,1)*GAy(5)*GAx(5))-(2*k5*sin(betta)*sin(ro)*f5(1,1)*f5(3,1)*GAy(5)*GAx(5))-(2*k5*sin(betta)*cos(ro)*f5(1,1)*f5(2,1)*GAy(5)*GAx(5))+(2*k5*cos(betta)*cos(ro)*f5(1,1)*f5(3,1)*GAy(5)*GAy(5))-(2*k5*sin(betta)*cos(ro)*f5(1,1)*f5(2,1)*GAy(5)*GAx(5))-(2*k5*cos(betta)*cos(ro)*f5(1,1)*f5(2,1)*GAy(5)*GAz(5))-(2*k5*sin(betta)*sin(ro)*f5(1,1)*f5(2,1)*GAx(5)*GAx(5))+(2*k5*cos(betta)*sin(ro)*f5(1,1)*f5(3,1)*GAy(5)*GAx(5))+(2*k5*sin(betta)*cos(ro)*f5(2,1)*f5(2,1)*GAx(5)*GAx(5))+(2*k5*sin(betta)*sin(ro)*f5(2,1)*f5(3,1)*GAx(5)*GAx(5))-(2*k5*cos(betta)*cos(ro)*f5(2,1)*f5(3,1)*GAx(5)*GAy(5))+(2*k5*cos(betta)*cos(ro)*f5(2,1)*f5(2,1)*GAz(5)*GAx(5))+(2*k5*cos(betta)*sin(ro)*f5(3,1)*f5(2,1)*GAz(5)*GAx(5))-(2*k5*cos(betta)*sin(ro)*f5(1,1)*f5(2,1)*GAz(5)*GAx(5))-(2*k5*cos(betta)*sin(ro)*f5(3,1)*f5(3,1)*GAy(5)*GAx(5))+(2*k5*sin(betta)*cos(ro)*f5(1,1)*f5(1,1)*GAy(5)*GAy(5))-(2*h5*sin(betta)*cos(ro)*GAx(5)*GAx(5))+(-2*h6*sin(betta)*cos(ro)*GAy(6)*GAy(6))-(2*h6*sin(betta)*sin(ro)*GAy(6)*GAz(6))-(2*h6*cos(betta)*cos(ro)*GAx(6)*GAz(6))+(2*h6*cos(betta)*sin(ro)*GAy(6)*GAx(6))+(2*k6*sin(betta)*sin(ro)*f6(1,1)*f6(1,1)*GAy(6)*GAx(6))-(2*k6*sin(betta)*sin(ro)*f6(1,1)*f6(3,1)*GAy(6)*GAx(6))-(2*k6*sin(betta)*cos(ro)*f6(1,1)*f6(2,1)*GAy(6)*GAx(6))+(2*k6*cos(betta)*cos(ro)*f6(1,1)*f6(3,1)*GAy(6)*GAy(6))-(2*k6*sin(betta)*cos(ro)*f6(1,1)*f6(2,1)*GAy(6)*GAx(6))-(2*k6*cos(betta)*cos(ro)*f6(1,1)*f6(2,1)*GAy(6)*GAz(6))-(2*k6*sin(betta)*sin(ro)*f6(1,1)*f6(2,1)*GAx(6)*GAx(6))+(2*k6*cos(betta)*sin(ro)*f6(1,1)*f6(3,1)*GAy(6)*GAx(6))+(2*k6*sin(betta)*cos(ro)*f6(2,1)*f6(2,1)*GAx(6)*GAx(6))+(2*k6*sin(betta)*sin(ro)*f6(2,1)*f6(3,1)*GAx(6)*GAx(6))-(2*k6*cos(betta)*cos(ro)*f6(2,1)*f6(3,1)*GAx(6)*GAy(6))+(2*k6*cos(betta)*cos(ro)*f6(2,1)*f6(2,1)*GAz(6)*GAx(6))+(2*k6*cos(betta)*sin(ro)*f6(3,1)*f6(2,1)*GAz(6)*GAx(6))-(2*k6*cos(betta)*sin(ro)*f6(1,1)*f6(2,1)*GAz(6)*GAx(6))-(2*k6*cos(betta)*sin(ro)*f6(3,1)*f6(3,1)*GAy(6)*GAx(6))+(2*k6*sin(betta)*cos(ro)*f6(1,1)*f6(1,1)*GAy(6)*GAy(6))-(2*h6*sin(betta)*cos(ro)*GAx(6)*GAx(6));
Mbegadot=(-2*h1*cos(ro)*GAy(1)*GAz(1))-(2*h1*sin(ro)*GAx(1)*GAx(1))-(2*h1*sin(ro)*GAz(1)*GAz(1))+(2*k1*cos(ro)*f1(1,1)*f1(1,1)*GAz(1)*GAy(1))+(2*k1*cos(ro)*f1(1,1)*f1(3,1)*GAx(1)*GAy(1))-(2*k1*cos(ro)*f1(2,1)*f1(1,1)*GAz(1)*GAx(1))-(2*k1*sin(ro)*f1(1,1)*f1(3,1)*GAz(1)*GAx(1))-(2*k1*sin(ro)*f1(1,1)*f1(3,1)*GAx(1)*GAx(1))+(2*k1*cos(ro)*f1(2,1)*f1(3,1)*GAx(1)*GAx(1))+(2*k1*sin(ro)*f1(3,1)*f1(3,1)*GAx(1)*GAx(1))+(2*k1*sin(ro)*f1(1,1)*f1(1,1)*GAz(1)*GAx(1))+(-2*h2*cos(ro)*GAy(2)*GAz(2))-(2*h2*sin(ro)*GAx(2)*GAx(2))-(2*h2*sin(ro)*GAz(2)*GAz(2))+(2*k2*cos(ro)*f2(1,1)*f2(1,1)*GAz(2)*GAy(2))+(2*k2*cos(ro)*f2(1,1)*f2(3,1)*GAx(2)*GAy(2))-(2*k2*cos(ro)*f2(2,1)*f2(1,1)*GAz(2)*GAx(2))-(2*k2*sin(ro)*f2(1,1)*f2(3,1)*GAz(2)*GAx(2))-(2*k2*sin(ro)*f2(1,1)*f2(3,1)*GAx(2)*GAx(2))+(2*k2*cos(ro)*f2(2,1)*f2(3,1)*GAx(2)*GAx(2))+(2*k2*sin(ro)*f2(3,1)*f2(3,1)*GAx(2)*GAx(2))+(2*k2*sin(ro)*f2(1,1)*f2(1,1)*GAz(2)*GAx(2))+(-2*h3*cos(ro)*GAy(3)*GAz(3))-(2*h3*sin(ro)*GAx(3)*GAx(3))-(2*h3*sin(ro)*GAz(3)*GAz(3))+(2*k3*cos(ro)*f3(1,1)*f3(1,1)*GAz(3)*GAy(3))+(2*k3*cos(ro)*f3(1,1)*f3(3,1)*GAx(3)*GAy(3))-(2*k3*cos(ro)*f3(2,1)*f3(1,1)*GAz(3)*GAx(3))-(2*k3*sin(ro)*f3(1,1)*f3(3,1)*GAz(3)*GAx(3))-(2*k3*sin(ro)*f3(1,1)*f3(3,1)*GAx(3)*GAx(3))+(2*k3*cos(ro)*f3(2,1)*f3(3,1)*GAx(3)*GAx(3))+(2*k3*sin(ro)*f3(3,1)*f3(3,1)*GAx(3)*GAx(3))+(2*k3*sin(ro)*f3(1,1)*f3(1,1)*GAz(3)*GAx(3))+(-2*h4*cos(ro)*GAy(4)*GAz(4))-(2*h4*sin(ro)*GAx(4)*GAx(4))-(2*h4*sin(ro)*GAz(4)*GAz(4))+(2*k4*cos(ro)*f4(1,1)*f4(1,1)*GAz(4)*GAy(4))+(2*k4*cos(ro)*f4(1,1)*f4(3,1)*GAx(4)*GAy(4))-(2*k4*cos(ro)*f4(2,1)*f4(1,1)*GAz(4)*GAx(4))-(2*k4*sin(ro)*f4(1,1)*f4(3,1)*GAz(4)*GAx(4))-(2*k4*sin(ro)*f4(1,1)*f4(3,1)*GAx(4)*GAx(4))+(2*k4*cos(ro)*f4(2,1)*f4(3,1)*GAx(4)*GAx(4))+(2*k4*sin(ro)*f4(3,1)*f4(3,1)*GAx(4)*GAx(4))+(2*k4*sin(ro)*f4(1,1)*f4(1,1)*GAz(4)*GAx(4))+(-2*h5*cos(ro)*GAy(5)*GAz(5))-(2*h5*sin(ro)*GAx(5)*GAx(5))-(2*h5*sin(ro)*GAz(5)*GAz(5))+(2*k5*cos(ro)*f5(1,1)*f5(1,1)*GAz(5)*GAy(5))+(2*k5*cos(ro)*f5(1,1)*f5(3,1)*GAx(5)*GAy(5))-(2*k5*cos(ro)*f5(2,1)*f5(1,1)*GAz(5)*GAx(5))-(2*k5*sin(ro)*f5(1,1)*f5(3,1)*GAz(5)*GAx(5))-(2*k5*sin(ro)*f5(1,1)*f5(3,1)*GAx(5)*GAx(5))+(2*k5*cos(ro)*f5(2,1)*f5(3,1)*GAx(5)*GAx(5))+(2*k5*sin(ro)*f5(3,1)*f5(3,1)*GAx(5)*GAx(5))+(2*k5*sin(ro)*f5(1,1)*f5(1,1)*GAz(5)*GAx(5))+(-2*h6*cos(ro)*GAy(6)*GAz(6))-(2*h6*sin(ro)*GAx(6)*GAx(6))-(2*h6*sin(ro)*GAz(6)*GAz(6))+(2*k6*cos(ro)*f6(1,1)*f6(1,1)*GAz(6)*GAy(6))+(2*k6*cos(ro)*f6(1,1)*f6(3,1)*GAx(6)*GAy(6))-(2*k6*cos(ro)*f6(2,1)*f6(1,1)*GAz(6)*GAx(6))-(2*k6*sin(ro)*f6(1,1)*f6(3,1)*GAz(6)*GAx(6))-(2*k6*sin(ro)*f6(1,1)*f6(3,1)*GAx(6)*GAx(6))+(2*k6*cos(ro)*f6(2,1)*f6(3,1)*GAx(6)*GAx(6))+(2*k6*sin(ro)*f6(3,1)*f6(3,1)*GAx(6)*GAx(6))+(2*k6*sin(ro)*f6(1,1)*f6(1,1)*GAz(6)*GAx(6));
MPxaldot=(2*h1*sin(betta)*GAy(1))+(2*k1*cos(betta)*f1(1,1)*f1(2,1)*GAz(1))+(2*k1*sin(betta)*f1(1,1)*f1(2,1)*GAx(1))-(2*k1*cos(betta)*f1(1,1)*f1(3,1)*GAy(1))-(2*k1*sin(betta)*f1(1,1)*f1(1,1)*GAy(1))+(2*h2*sin(betta)*GAy(2))+(2*k2*cos(betta)*f2(1,1)*f2(2,1)*GAz(2))+(2*k2*sin(betta)*f2(1,1)*f2(2,1)*GAx(2))-(2*k2*cos(betta)*f2(1,1)*f2(3,1)*GAy(2))-(2*k2*sin(betta)*f2(1,1)*f2(1,1)*GAy(2))+(2*h3*sin(betta)*GAy(3))+(2*k3*cos(betta)*f3(1,1)*f3(2,1)*GAz(3))+(2*k3*sin(betta)*f3(1,1)*f3(2,1)*GAx(3))-(2*k3*cos(betta)*f3(1,1)*f3(3,1)*GAy(3))-(2*k3*sin(betta)*f3(1,1)*f3(1,1)*GAy(3))+(2*h4*sin(betta)*GAy(4))+(2*k4*cos(betta)*f4(1,1)*f4(2,1)*GAz(4))+(2*k4*sin(betta)*f4(1,1)*f4(2,1)*GAx(4))-(2*k4*cos(betta)*f4(1,1)*f4(3,1)*GAy(4))-(2*k4*sin(betta)*f4(1,1)*f4(1,1)*GAy(4))+(2*h5*sin(betta)*GAy(5))+(2*k5*cos(betta)*f5(1,1)*f5(2,1)*GAz(5))+(2*k5*sin(betta)*f5(1,1)*f5(2,1)*GAx(5))-(2*k5*cos(betta)*f5(1,1)*f5(3,1)*GAy(5))-(2*k5*sin(betta)*f5(1,1)*f5(1,1)*GAy(5))+(2*h6*sin(betta)*GAy(6))+(2*k6*cos(betta)*f6(1,1)*f6(2,1)*GAz(6))+(2*k6*sin(betta)*f6(1,1)*f6(2,1)*GAx(6))-(2*k6*cos(betta)*f6(1,1)*f6(3,1)*GAy(6))-(2*k6*sin(betta)*f6(1,1)*f6(1,1)*GAy(6));
MPxbedot=(2*h1*GAz(1))-(2*k1*f1(1,1)*f1(1,1)*GAz(1))+(2*k1*f1(1,1)*f1(3,1)*GAx(1))+(2*h2*GAz(2))-(2*k2*f2(1,1)*f2(1,1)*GAz(2))+(2*k2*f2(1,1)*f2(3,1)*GAx(2))+(2*h3*GAz(3))-(2*k3*f3(1,1)*f3(1,1)*GAz(3))+(2*k3*f3(1,1)*f3(3,1)*GAx(3))+(2*h4*GAz(4))-(2*k4*f4(1,1)*f4(1,1)*GAz(4))+(2*k4*f4(1,1)*f4(3,1)*GAx(4))+(2*h5*GAz(5))-(2*k5*f5(1,1)*f5(1,1)*GAz(5))+(2*k5*f5(1,1)*f5(3,1)*GAx(5))+(2*h6*GAz(6))-(2*k6*f6(1,1)*f6(1,1)*GAz(6))+(2*k6*f6(1,1)*f6(3,1)*GAx(6));
MPxgadot=(2*h1*cos(ro)*GAy(1))-(2*h1*sin(ro)*GAz(1))-(2*k1*cos(ro)*f1(1,1)*f1(2,1)*GAx(1))+(2*k1*sin(ro)*f1(1,1)*f1(1,1)*GAx(1))+(2*k1*cos(ro)*f1(1,1)*f1(1,1)*GAy(1))-(2*k1*sin(ro)*f1(1,1)*f1(3,1)*GAx(1))+(2*h2*cos(ro)*GAy(2))-(2*h2*sin(ro)*GAz(2))-(2*k2*cos(ro)*f2(1,1)*f2(2,1)*GAx(2))+(2*k2*sin(ro)*f2(1,1)*f2(1,1)*GAx(2))+(2*k2*cos(ro)*f2(1,1)*f2(1,1)*GAy(2))-(2*k2*sin(ro)*f2(1,1)*f2(3,1)*GAx(2))+(2*h3*cos(ro)*GAy(3))-(2*h3*sin(ro)*GAz(3))-(2*k3*cos(ro)*f3(1,1)*f3(2,1)*GAx(3))+(2*k3*sin(ro)*f3(1,1)*f3(1,1)*GAx(3))+(2*k3*cos(ro)*f3(1,1)*f3(1,1)*GAy(3))-(2*k3*sin(ro)*f3(1,1)*f3(3,1)*GAx(3))+(2*h4*cos(ro)*GAy(4))-(2*h4*sin(ro)*GAz(4))-(2*k4*cos(ro)*f4(1,1)*f4(2,1)*GAx(4))+(2*k4*sin(ro)*f4(1,1)*f4(1,1)*GAx(4))+(2*k4*cos(ro)*f4(1,1)*f4(1,1)*GAy(4))-(2*k4*sin(ro)*f4(1,1)*f4(3,1)*GAx(4))+(2*h5*cos(ro)*GAy(5))-(2*h5*sin(ro)*GAz(5))-(2*k5*cos(ro)*f5(1,1)*f5(2,1)*GAx(5))+(2*k5*sin(ro)*f5(1,1)*f5(1,1)*GAx(5))+(2*k5*cos(ro)*f5(1,1)*f5(1,1)*GAy(5))-(2*k5*sin(ro)*f5(1,1)*f5(3,1)*GAx(5))+(2*h6*cos(ro)*GAy(6))-(2*h6*sin(ro)*GAz(6))-(2*k6*cos(ro)*f6(1,1)*f6(2,1)*GAx(6))+(2*k6*sin(ro)*f6(1,1)*f6(1,1)*GAx(6))+(2*k6*cos(ro)*f6(1,1)*f6(1,1)*GAy(6))-(2*k6*sin(ro)*f6(1,1)*f6(3,1)*GAx(6));
MPxPydot=(-2*k1*f1(1,1)*f1(2,1))+(-2*k2*f2(1,1)*f2(2,1))+(-2*k3*f3(1,1)*f3(2,1))+(-2*k4*f4(1,1)*f4(2,1))+(-2*k5*f5(1,1)*f5(2,1))+(-2*k6*f6(1,1)*f6(2,1));
MPxPzdot=(-2*k1*f1(1,1)*f1(3,1))+(-2*k2*f2(1,1)*f2(3,1))+(-2*k3*f3(1,1)*f3(3,1))+(-2*k4*f4(1,1)*f4(3,1))+(-2*k5*f5(1,1)*f5(3,1))+(-2*k6*f6(1,1)*f6(3,1));
MPyaldot=(-2*h1*cos(betta)*GAz(1))-(2*h1*sin(betta)*GAx(1))-(2*k1*sin(betta)*f1(1,1)*f1(2,1)*GAy(1))+(2*k1*sin(betta)*f1(2,1)*f1(2,1)*GAx(1))+(2*k1*cos(betta)*f1(2,1)*f1(2,1)*GAz(1))-(2*k1*cos(betta)*f1(2,1)*f1(3,1)*GAy(1))+(-2*h2*cos(betta)*GAz(2))-(2*h2*sin(betta)*GAx(2))-(2*k2*sin(betta)*f2(1,1)*f2(2,1)*GAy(2))+(2*k2*sin(betta)*f2(2,1)*f2(2,1)*GAx(2))+(2*k2*cos(betta)*f2(2,1)*f2(2,1)*GAz(2))-(2*k2*cos(betta)*f2(2,1)*f2(3,1)*GAy(2))+(-2*h3*cos(betta)*GAz(3))-(2*h3*sin(betta)*GAx(3))-(2*k3*sin(betta)*f3(1,1)*f3(2,1)*GAy(3))+(2*k3*sin(betta)*f3(2,1)*f3(2,1)*GAx(3))+(2*k3*cos(betta)*f3(2,1)*f3(2,1)*GAz(3))-(2*k3*cos(betta)*f3(2,1)*f3(3,1)*GAy(3))+(-2*h4*cos(betta)*GAz(4))-(2*h4*sin(betta)*GAx(4))-(2*k4*sin(betta)*f4(1,1)*f4(2,1)*GAy(4))+(2*k4*sin(betta)*f4(2,1)*f4(2,1)*GAx(4))+(2*k4*cos(betta)*f4(2,1)*f4(2,1)*GAz(4))-(2*k4*cos(betta)*f4(2,1)*f4(3,1)*GAy(4))+(-2*h5*cos(betta)*GAz(5))-(2*h5*sin(betta)*GAx(5))-(2*k5*sin(betta)*f5(1,1)*f5(2,1)*GAy(5))+(2*k5*sin(betta)*f5(2,1)*f5(2,1)*GAx(5))+(2*k5*cos(betta)*f5(2,1)*f5(2,1)*GAz(5))-(2*k5*cos(betta)*f5(2,1)*f5(3,1)*GAy(5))+(-2*h6*cos(betta)*GAz(6))-(2*h6*sin(betta)*GAx(6))-(2*k6*sin(betta)*f6(1,1)*f6(2,1)*GAy(6))+(2*k6*sin(betta)*f6(2,1)*f6(2,1)*GAx(6))+(2*k6*cos(betta)*f6(2,1)*f6(2,1)*GAz(6))-(2*k6*cos(betta)*f6(2,1)*f6(3,1)*GAy(6));
MPybedot=(-2*k1*f1(1,1)*f1(2,1)*GAz(1))+(2*k1*f1(3,1)*f1(2,1)*GAx(1))+(-2*k1*f2(1,1)*f2(2,1)*GAz(2))+(2*k2*f2(3,1)*f2(2,1)*GAx(2))+(-2*k1*f3(1,1)*f3(2,1)*GAz(3))+(2*k3*f3(3,1)*f3(2,1)*GAx(3))+(-2*k1*f4(1,1)*f4(2,1)*GAz(4))+(2*k4*f4(3,1)*f4(2,1)*GAx(4))+(-2*k1*f5(1,1)*f5(2,1)*GAz(5))+(2*k5*f5(3,1)*f5(2,1)*GAx(5))+(-2*k1*f6(1,1)*f6(2,1)*GAz(6))+(2*k6*f6(3,1)*f6(2,1)*GAx(6));
MPygadot=(2*h1*cos(ro)*GAx(1))+(2*k1*cos(ro)*f1(2,1)*f1(1,1)*GAy(1))+(2*k1*sin(ro)*f1(2,1)*f1(1,1)*GAx(1))-(2*k1*cos(ro)*f1(2,1)*f1(2,1)*GAx(1))-(2*k1*sin(ro)*f1(2,1)*f1(3,1)*GAx(1))+(2*h2*cos(ro)*GAx(2))+(2*k2*cos(ro)*f2(2,1)*f2(1,1)*GAy(2))+(2*k2*sin(ro)*f2(2,1)*f2(1,1)*GAx(2))-(2*k2*cos(ro)*f2(2,1)*f2(2,1)*GAx(2))-(2*k2*sin(ro)*f2(2,1)*f2(3,1)*GAx(2))+(2*h3*cos(ro)*GAx(3))+(2*k3*cos(ro)*f3(2,1)*f3(1,1)*GAy(3))+(2*k3*sin(ro)*f3(2,1)*f3(1,1)*GAx(3))-(2*k3*cos(ro)*f3(2,1)*f3(2,1)*GAx(3))-(2*k3*sin(ro)*f3(2,1)*f3(3,1)*GAx(3))+(2*h4*cos(ro)*GAx(4))+(2*k4*cos(ro)*f4(2,1)*f4(1,1)*GAy(4))+(2*k4*sin(ro)*f4(2,1)*f4(1,1)*GAx(4))-(2*k4*cos(ro)*f4(2,1)*f4(2,1)*GAx(4))-(2*k4*sin(ro)*f4(2,1)*f4(3,1)*GAx(4))+(2*h5*cos(ro)*GAx(5))+(2*k5*cos(ro)*f5(2,1)*f5(1,1)*GAy(5))+(2*k5*sin(ro)*f5(2,1)*f5(1,1)*GAx(5))-(2*k5*cos(ro)*f5(2,1)*f5(2,1)*GAx(5))-(2*k5*sin(ro)*f5(2,1)*f5(3,1)*GAx(5))+(2*h6*cos(ro)*GAx(6))+(2*k6*cos(ro)*f6(2,1)*f6(1,1)*GAy(6))+(2*k6*sin(ro)*f6(2,1)*f6(1,1)*GAx(6))-(2*k6*cos(ro)*f6(2,1)*f6(2,1)*GAx(6))-(2*k6*sin(ro)*f6(2,1)*f6(3,1)*GAx(6));
MPyPzdot=(-2*k1*f1(2,1)*f1(3,1))+(-2*k2*f2(2,1)*f2(3,1))+(-2*k3*f3(2,1)*f3(3,1))+(-2*k4*f4(2,1)*f4(3,1))+(-2*k5*f5(2,1)*f5(3,1))+(-2*k6*f6(2,1)*f6(3,1));
MPzaldot=(-2*k1*sin(betta)*f1(3,1)*f1(1,1)*GAy(1))+(2*k1*sin(betta)*f1(2,1)*f1(3,1)*GAx(1))+(2*k1*cos(betta)*f1(2,1)*f1(3,1)*GAz(1))-(2*k1*cos(betta)*f1(3,1)*f1(3,1)*GAy(1))+(2*h1*cos(betta)*GAy(1))+(-2*k2*sin(betta)*f2(3,1)*f2(1,1)*GAy(2))+(2*k2*sin(betta)*f2(2,1)*f2(3,1)*GAx(2))+(2*k2*cos(betta)*f2(2,1)*f2(3,1)*GAz(2))-(2*k2*cos(betta)*f2(3,1)*f2(3,1)*GAy(2))+(2*h2*cos(betta)*GAy(2))+(-2*k3*sin(betta)*f3(3,1)*f3(1,1)*GAy(3))+(2*k3*sin(betta)*f3(2,1)*f3(3,1)*GAx(3))+(2*k3*cos(betta)*f3(2,1)*f3(3,1)*GAz(3))-(2*k3*cos(betta)*f3(3,1)*f3(3,1)*GAy(3))+(2*h3*cos(betta)*GAy(3))+(-2*k4*sin(betta)*f4(3,1)*f4(1,1)*GAy(4))+(2*k4*sin(betta)*f4(2,1)*f4(3,1)*GAx(4))+(2*k4*cos(betta)*f4(2,1)*f4(3,1)*GAz(4))-(2*k4*cos(betta)*f4(3,1)*f4(3,1)*GAy(4))+(2*h4*cos(betta)*GAy(4))+(-2*k5*sin(betta)*f5(3,1)*f5(1,1)*GAy(5))+(2*k5*sin(betta)*f5(2,1)*f5(3,1)*GAx(5))+(2*k5*cos(betta)*f5(2,1)*f5(3,1)*GAz(5))-(2*k5*cos(betta)*f5(3,1)*f5(3,1)*GAy(5))+(2*h5*cos(betta)*GAy(5))+(-2*k6*sin(betta)*f6(3,1)*f6(1,1)*GAy(6))+(2*k6*sin(betta)*f6(2,1)*f6(3,1)*GAx(6))+(2*k6*cos(betta)*f6(2,1)*f6(3,1)*GAz(6))-(2*k6*cos(betta)*f6(3,1)*f6(3,1)*GAy(6))+(2*h6*cos(betta)*GAy(6));
MPzbedot=(-2*h1*GAx(1))-(2*k1*f1(3,1)*f1(1,1)*GAx(1))+(2*k1*f1(3,1)*f1(3,1)*GAx(1))+(-2*h2*GAx(2))-(2*k2*f2(3,1)*f2(1,1)*GAx(2))+(2*k2*f1(3,1)*f2(3,1)*GAx(2))+(-2*h3*GAx(3))-(2*k3*f3(3,1)*f3(1,1)*GAx(3))+(2*k3*f1(3,1)*f3(3,1)*GAx(3))+(-2*h4*GAx(4))-(2*k4*f4(3,1)*f4(1,1)*GAx(4))+(2*k4*f1(3,1)*f4(3,1)*GAx(4))+(-2*h5*GAx(5))-(2*k5*f5(3,1)*f5(1,1)*GAx(5))+(2*k5*f1(3,1)*f5(3,1)*GAx(5))+(-2*h6*GAx(6))-(2*k6*f6(3,1)*f6(1,1)*GAx(6))+(2*k6*f1(3,1)*f6(3,1)*GAx(6));
MPzgadot=(-2*k1*sin(ro)*f1(3,1)*f1(3,1)*GAx(1))-(2*k1*cos(ro)*f1(3,1)*f1(2,1)*GAx(1))+(2*k1*cos(ro)*f1(3,1)*f1(1,1)*GAy(1))+(2*k1*sin(ro)*f1(3,1)*f1(1,1)*GAx(1))+(2*h1*sin(ro)*GAx(1))+(-2*k2*sin(ro)*f2(3,1)*f2(3,1)*GAx(2))-(2*k2*cos(ro)*f2(3,1)*f2(2,1)*GAx(2))+(2*k2*cos(ro)*f2(3,1)*f2(1,1)*GAy(2))+(2*k2*sin(ro)*f2(3,1)*f2(1,1)*GAx(2))+(2*h2*sin(ro)*GAx(2))+(-2*k3*sin(ro)*f3(3,1)*f3(3,1)*GAx(3))-(2*k3*cos(ro)*f3(3,1)*f3(2,1)*GAx(3))+(2*k3*cos(ro)*f3(3,1)*f3(1,1)*GAy(3))+(2*k3*sin(ro)*f3(3,1)*f3(1,1)*GAx(3))+(2*h3*sin(ro)*GAx(3))+(-2*k4*sin(ro)*f4(3,1)*f4(3,1)*GAx(4))-(2*k4*cos(ro)*f4(3,1)*f4(2,1)*GAx(4))+(2*k4*cos(ro)*f4(3,1)*f4(1,1)*GAy(4))+(2*k4*sin(ro)*f4(3,1)*f4(1,1)*GAx(4))+(2*h4*sin(ro)*GAx(4))+(-2*k5*sin(ro)*f5(3,1)*f5(3,1)*GAx(5))-(2*k5*cos(ro)*f5(3,1)*f5(2,1)*GAx(5))+(2*k5*cos(ro)*f5(3,1)*f5(1,1)*GAy(5))+(2*k5*sin(ro)*f5(3,1)*f5(1,1)*GAx(5))+(2*h5*sin(ro)*GAx(5))+(-2*k6*sin(ro)*f6(3,1)*f6(3,1)*GAx(6))-(2*k6*cos(ro)*f6(3,1)*f6(2,1)*GAx(6))+(2*k6*cos(ro)*f6(3,1)*f6(1,1)*GAy(6))+(2*k6*sin(ro)*f6(3,1)*f6(1,1)*GAx(6))+(2*h6*sin(ro)*GAx(6));
Mlegs=.5*(m1+m2)*[MPxdot ,.5*MPxPydot,.5*MPxPzdot,.5*MPxaldot,.5*MPxbedot,.5*MPxgadot;.5*MPxPydot,MPydot ,.5*MPyPzdot,.5*MPyaldot,.5*MPybedot,.5*MPygadot;.5*MPxPzdot,.5*MPyPzdot,MPzdot ,.5*MPzaldot,.5*MPzbedot,.5*MPzgadot;.5*MPxaldot,.5*MPyaldot,.5*MPzaldot,Maldot ,.5*Malbedot,.5*Malgadot;.5*MPxbedot,.5*MPybedot,.5*MPzbedot,.5*Malbedot,Mbedot ,.5*Mbegadot;.5*MPxgadot,.5*MPygadot,.5*MPzgadot,.5*Malgadot,.5*Mbegadot,Mgadot ]
disp('Mlegs:'),disp(Mlegs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Kinetic energy: legs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Potential energy: legs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1
ZA1=[0,0,1]*Rz*Rx*Ry*[GAx(j);GAy(j);GAz(j)];
disp('ZA1:'),disp(ZA1)
end
for j=2
ZA2=[0,0,1]*Rz*Rx*Ry*[GAx(j);GAy(j);GAz(j)];
disp('ZA2:'),disp(ZA2)
end
for j=3
ZA3=[0,0,1]*Rz*Rx*Ry*[GAx(j);GAy(j);GAz(j)];
disp('ZA3:'),disp(ZA3)
end
Plegs=(m1+m2)*g*(((Ia * ((1/l2)+(1/l1))) + ((2*m2)/(m1+m2)) * (Pz+ZA1)) + ((Ia * ((1/l4)+(1/l3))) + ((2*m2)/(m1+m2) ) * (Pz+ZA2)) + ((Ia * ((1/l6)+(1/l5))) + ((2*m2)/(m1+m2) ) * (Pz+ZA3)) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Potential energy: legs %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START Dynamic equations %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START Vm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MXp = Mpl+Mlegs;
disp('MXp:'),disp(MXp);
for i=1:6
for j=1:6
Vmpl=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=1
for j=1
Vmpl11=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) +((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vmpl12=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vmpl13=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vmpl14=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vmpl15=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vmpl16=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=2
for j=1
Vmpl21=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vmpl22=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) +((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vmpl23=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vmpl24=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vmpl25=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) +((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vmpl26=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=3
for j=1
Vmpl31=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vmpl32=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3;Vmup33=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vmpl34=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vmpl35=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vmpl36=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=4
for j=1
Vmpl41=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vmpl42=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vmpl43=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vmpl44=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vmpl45=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vmpl46=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=5
for j=1
Vmpl51=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) +((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vmpl52=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vmpl53=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vmpl54=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) +((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vmpl55=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vmpl56=0.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=6
for j=1
Vmpl61=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vmpl62=0.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vmpl63=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vmpl64=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vmpl65=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vmpl66=.5*((diff(Mpl(1,j),Xp(i,1)))+(diff(Mpl(1,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mpl(2,j),Xp(i,1)))+(diff(Mpl(2,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mpl(3,j),Xp(i,1)))+(diff(Mpl(3,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mpl(4,j),Xp(i,1)))+(diff(Mpl(4,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mpl(5,j),Xp(i,1)))+(diff(Mpl(5,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mpl(6,j),Xp(i,1)))+(diff(Mpl(6,i),Xp(j,1)))+(diff(Mpl(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
Vmpl=[Vmpl11,Vmpl12,Vmpl13,Vmpl14,Vmpl15,Vmpl16;
Vmpl21,Vmpl22,Vmpl23,Vmpl24,Vmpl25,Vmpl26;
Vmpl31,Vmpl32,Vmup33,Vmpl34,Vmpl35,Vmpl36;
Vmpl41,Vmpl42,Vmpl43,Vmpl44,Vmpl45,Vmpl46;
Vmpl51,Vmpl52,Vmpl53,Vmpl54,Vmpl55,Vmpl56;
Vmpl61,Vmpl62,Vmpl63,Vmpl64,Vmpl65,Vmpl66]
for i=1
for j=1
Vlegs11=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vlegs12=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vlegs13=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vlegs14=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vlegs15=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) +((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vlegs16=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=2
for j=1
Vlegs21=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vlegs22=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vlegs23=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vlegs24=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vlegs25=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
  Vlegs26=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=3
for j=1
Vlegs31=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vlegs32=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vlegs33=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) +  ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vlegs34=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vlegs35=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vlegs36=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=4
for j=1
Vlegs41=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vlegs42=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vlegs43=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vlegs44=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) +((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vlegs45=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vlegs46=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=5
for j=1
Vlegs51=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vlegs52=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vlegs53=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vlegs54=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
    Vlegs55=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vlegs56=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
for i=6
for j=1
Vlegs61=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=2
Vlegs62=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) +((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=3
Vlegs63=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=4
Vlegs64=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=5
Vlegs65=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
for j=6
Vlegs66=.5*((diff(Mlegs(1,j),Xp(i,1)))+(diff(Mlegs(1,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(1,1)))*(Xpdot(1,1))) + ((diff(Mlegs(2,j),Xp(i,1)))+(diff(Mlegs(2,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(2,1)))*(Xpdot(2,1))) + ((diff(Mlegs(3,j),Xp(i,1)))+(diff(Mlegs(3,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(3,1)))*(Xpdot(3,1))) + ((diff(Mlegs(4,j),Xp(i,1)))+(diff(Mlegs(4,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(4,1)))*(Xpdot(4,1))) + ((diff(Mlegs(5,j),Xp(i,1)))+(diff(Mlegs(5,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(5,1)))*(Xpdot(5,1))) + ((diff(Mlegs(6,j),Xp(i,1)))+(diff(Mlegs(6,i),Xp(j,1)))+(diff(Mlegs(i,j),Xp(6,1)))*(Xpdot(6,1)));
end
end
Vlegs=[Vlegs11,Vlegs12,Vlegs13,Vlegs14,Vlegs15,Vlegs16;
Vlegs21,Vlegs22,Vlegs23,Vlegs24,Vlegs25,Vlegs26;
Vlegs31,Vlegs32,Vlegs33,Vlegs34,Vlegs35,Vlegs36;
Vlegs41,Vlegs42,Vlegs43,Vlegs44,Vlegs45,Vlegs46;
Vlegs51,Vlegs52,Vlegs53,Vlegs54,Vlegs55,Vlegs56;
Vlegs61,Vlegs62,Vlegs63,Vlegs64,Vlegs65,Vlegs66]
VXp = Vmpl+Vlegs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Vm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START Gup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gpl11=diff((Ppl),Xp(1,1));
Gpl21=diff((Ppl),Xp(2,1));
Gpl31=diff((Ppl),Xp(3,1));
Gpl41=diff((Ppl),Xp(4,1));
Gpl51=diff((Ppl),Xp(5,1));
Gpl61=diff((Ppl),Xp(6,1));
Gpl = [Gpl11;Gpl21;Gpl31;Gpl41;Gpl51;Gpl61]
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Gup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%START Glegs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Glegs11=diff((Plegs),Xp(1,1));
Glegs21=diff((Plegs),Xp(2,1));
Glegs31=diff((Plegs),Xp(3,1));
Glegs41=diff((Plegs),Xp(4,1));
Glegs51=diff((Plegs),Xp(5,1));
Glegs61=diff((Plegs),Xp(6,1));
Glegs = [Glegs11;Glegs21;Glegs31;Glegs41;Glegs51;Glegs61]
GXp = Gpl + Glegs
%%%%%%%%%%%%%%%%%%% END Glegs %%%%%%%%%%%%%%%
%syms Xpddot alddot beddot gaddot Xddot Yddot Zddot F f1 f2 f3 f4 f5 f6 jjjj
Xpddot = [Xddot,Yddot,Zddot,alddot,beddot,gaddot]'
JBB = JB'
% F = solve([JBB*F == ((MXp*Xpddot)+(VXp*Xpdot)+GXp)],[F])
% solve('JBB*F==((MXp*Xpddot)+(VXp*Xpdot)+GXp)','F')
F = JBB \ (((MXp*Xpddot)+(VXp*Xpdot)+GXp)); 
% Xpddot:6*1 MXp:6*6 JBB:6*6 VXp:6*6 Xpdot:6*1 GXp:6*1
disp('F:'),disp(F);
% ff1 = F(1,1)
% ff2 = F(2,1)
% ff3 = F(3,1)
% ff4 = F(4,1)
% ff5 = F(5,1)
% ff6 = F(6,1)

end
