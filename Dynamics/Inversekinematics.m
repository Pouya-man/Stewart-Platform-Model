function [l] = Inversekinematics( thetap,thetab,rp,rb,Px, PPy, Pz, ro, betta, fii, Pxdot, Pydot, Pzdot, rodot, bettadot, fiidot, mup, Ix, Iy, Iz, g, del, m1, m2, lj1, lj2 )
%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%
clc
clear all
syms thetap thetab rp rb ro betta fii Px PPy Pz rodot bettadot fiidot Pxdot Pydot Pzdot aldot bedot gadot mup Ix Iy Iz g del m1 m2 lj1 lj2 k GAxx GAyy GAzz wx wy wz ww Pdot tjj hi vv uix uiy uiz kii
P=[Px;PPy;Pz];
Pdot=[Pxdot;Pydot;Pzdot];
%%%%%%%%%%%%%%%%% Inverse kinematics %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% The connection points to the fixed platform (GA(i)) %%%%%%%%%%%%%%%%%
for i=[1,3,5]
landa(i)=((i*(pi/3))+(thetap/2)-(pi/3));
disp('landa:'), disp(i),disp(landa(i));
end
for i=[2,4,6]
landa(i)=landa(i-1)-thetap+(((2*pi)/3));
disp('landa:'), disp(i),disp(landa(i))
end
for i=1:6
GAx(i)=rp*cos(landa(i));
GAy(i)=rp*sin(landa(i));
GAz(i)=0;
end
for i=1:6
GA=[GAx(i);GAy(i);GAz(i)];
disp('GA:'),disp(i),disp(GA)
end
%%%%%%%%%%%%%%%%% The connection points to the moving platform (B(i)) %%%%%%%%%%%%%%%%%
for i=[1,3,5]
delta(i)=((i*(pi/3))-(thetab/2));
disp('delta:'), disp(i),disp(delta(i))
end
for i=[2,4,6]
delta(i)=delta(i-1)+ thetab;
disp('delta:'), disp(i),disp(delta(i))
end
for i=1:6
Bx(i)=rb*cos(delta(i));
By(i)=rb*sin(delta(i));
Bz(i)=0;
end
for i=1:6
B=[Bx(i);By(i);Bz(i)];
disp('B:'),disp(i),disp(B)
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
disp('BRA:'),disp(BRA)
%%%%%%%%%%%%%%%%% END BRA %%%%%%%%%%%%%%%%%
ww=[0,-((cos(ro)*fiidot)+(-sin(betta)*rodot)),(bettadot-(fiidot*sin(ro)));((cos(ro)*fiidot)+(-sin(betta)*rodot)),0,-(cos(betta)*rodot);-(bettadot-(fiidot*sin(ro))),(cos(betta)*rodot),0]
T=(ww*([GAxx;GAyy;GAzz]))
Tjj=Pdot+T
vv=Tjj'*hi*Tjj
uii=[uix;uiy;uiz]
Vii=(Pdot'*uii)+(T'*uii)
cc=Vii*kii*Vii
ll=.5*(m1+m2)*(vv-cc)
%%%%%%%%%%%%%%%%% Li: link vector %%%%%%%%%%%%%%%%%
for i=1:6
L=(BRA*[GAx(i);GAy(i);GAz(i)])+P-[Bx(i);By(i);Bz(i)];
disp('L:'),disp(i),disp(L)
end
%%%%%%%%%%%%%%%%% END Li %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  l: length of the arms %%%%%%%%%%%%%%%%%
for i=1:6
l=(((Px-Bx(i))+(GAx(i)*r11)+(GAy(i)*r12))^2+((PPy-By(i))+(GAx(i)*r21)+(GAy(i)*r22))^2+((Pz)+(GAx(i)*r31)+(GAy(i)*r32))^2)^.5;
disp('l:'),disp(i),disp(l)
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