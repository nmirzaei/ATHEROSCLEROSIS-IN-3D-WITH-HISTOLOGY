%read the data Lumen1 and Intmed1
n = 250   %should be the same as the mesh file n's
x1 = 0;
x2 = 2*pi;

[a0,a,b] = seg(LumenA(:,1),LumenA(:,2),1);
[aa0,aa,bb] = seg(IntMedA(:,1),IntMedA(:,2),1);

sprintf('The lumen equation is (%f)+(%f)*Cos(x)+(%f)*Cos(2*x)+(%f)*Cos(3*x)+(%f)*Cos(4*x)+(%f)*Cos(5*x)+(%f)*Cos(6*x)+(%f)*Cos(7*x)+(%f)*Cos(8*x)+(%f)*Sin(x)+(%f)*Sin(2*x)+(%f)*Sin(3*x)+(%f)*Sin(4*x)+(%f)*Sin(5*x)+(%f)*Sin(6*x)+(%f)*Sin(7*x)+(%f)*Sin(8*x)',a0,a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8))

sprintf('The intmed equation is (%f)+(%f)*Cos(x)+(%f)*Cos(2*x)+(%f)*Cos(3*x)+(%f)*Cos(4*x)+(%f)*Cos(5*x)+(%f)*Cos(6*x)+(%f)*Cos(7*x)+(%f)*Cos(8*x)+(%f)*Sin(x)+(%f)*Sin(2*x)+(%f)*Sin(3*x)+(%f)*Sin(4*x)+(%f)*Sin(5*x)+(%f)*Sin(6*x)+(%f)*Sin(7*x)+(%f)*Sin(8*x)',aa0,aa(1),aa(2),aa(3),aa(4),aa(5),aa(6),aa(7),aa(8),bb(1),bb(2),bb(3),bb(4),bb(5),bb(6),bb(7),bb(8))


[A0,A,B] = seg(LumenB(:,1),LumenB(:,2),1);
[AA0,AA,BB] = seg(IntMedB(:,1),IntMedB(:,2),1);

sprintf('The lumen equation is (%f)+(%f)*Cos(x)+(%f)*Cos(2*x)+(%f)*Cos(3*x)+(%f)*Cos(4*x)+(%f)*Cos(5*x)+(%f)*Cos(6*x)+(%f)*Cos(7*x)+(%f)*Cos(8*x)+(%f)*Sin(x)+(%f)*Sin(2*x)+(%f)*Sin(3*x)+(%f)*Sin(4*x)+(%f)*Sin(5*x)+(%f)*Sin(6*x)+(%f)*Sin(7*x)+(%f)*Sin(8*x)',A0,A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),B(1),B(2),B(3),B(4),B(5),B(6),B(7),B(8))

sprintf('The intmed equation is (%f)+(%f)*Cos(x)+(%f)*Cos(2*x)+(%f)*Cos(3*x)+(%f)*Cos(4*x)+(%f)*Cos(5*x)+(%f)*Cos(6*x)+(%f)*Cos(7*x)+(%f)*Cos(8*x)+(%f)*Sin(x)+(%f)*Sin(2*x)+(%f)*Sin(3*x)+(%f)*Sin(4*x)+(%f)*Sin(5*x)+(%f)*Sin(6*x)+(%f)*Sin(7*x)+(%f)*Sin(8*x)',AA0,AA(1),AA(2),AA(3),AA(4),AA(5),AA(6),AA(7),AA(8),BB(1),BB(2),BB(3),BB(4),BB(5),BB(6),BB(7),BB(8))


% fA_int = @(x) (a0)+(a(1))*cos(x)+(a(2))*cos(2*x)+(a(3))*cos(3*x)+(a(4))*cos(4*x)+(a(5))*cos(5*x)+(a(6))*cos(6*x)+(a(7))*cos(7*x)+(a(8))*cos(8*x)+(b(1))*sin(x)+(b(2))*sin(2*x)+(b(3))*sin(3*x)+(b(4))*sin(4*x)+(b(5))*sin(5*x)+(b(6))*sin(6*x)+(b(7))*sin(7*x)+(b(8))*sin(8*x)
% fA_med = @(x) (aa0)+(aa(1))*cos(x)+(aa(2))*cos(2*x)+(aa(3))*cos(3*x)+(aa(4))*cos(4*x)+(aa(5))*cos(5*x)+(aa(6))*cos(6*x)+(aa(7))*cos(7*x)+(aa(8))*cos(8*x)+(bb(1))*sin(x)+(bb(2))*sin(2*x)+(bb(3))*sin(3*x)+(bb(4))*sin(4*x)+(bb(5))*sin(5*x)+(bb(6))*sin(6*x)+(bb(7))*sin(7*x)+(bb(8))*sin(8*x)
% 
% 
% fB_int = @(x) (A0)+(A(1))*cos(x)+(A(2))*cos(2*x)+(A(3))*cos(3*x)+(A(4))*cos(4*x)+(A(5))*cos(5*x)+(A(6))*cos(6*x)+(A(7))*cos(7*x)+(A(8))*cos(8*x)+(B(1))*sin(x)+(B(2))*sin(2*x)+(B(3))*sin(3*x)+(B(4))*sin(4*x)+(B(5))*sin(5*x)+(B(6))*sin(6*x)+(B(7))*sin(7*x)+(B(8))*sin(8*x)
% fB_med = @(x) (AA0)+(AA(1))*cos(x)+(AA(2))*cos(2*x)+(AA(3))*cos(3*x)+(AA(4))*cos(4*x)+(AA(5))*cos(5*x)+(AA(6))*cos(6*x)+(AA(7))*cos(7*x)+(AA(8))*cos(8*x)+(BB(1))*sin(x)+(BB(2))*sin(2*x)+(BB(3))*sin(3*x)+(BB(4))*sin(4*x)+(BB(5))*sin(5*x)+(BB(6))*sin(6*x)+(BB(7))*sin(7*x)+(BB(8))*sin(8*x)
% 
% for i= 0:n
% 
%   x = x1 + (x2 - x1) * i/n;
% 
%   r1_A(i+1) = fA_int(x); 
%   r2_A(i+1) = fA_med(x);
%   r1_B(i+1) = fB_int(x);
%   r2_B(i+1) = fB_med(x);
%   PointA_int(i+1,:) = [r1_A(i+1)*cos(x), r1_A(i+1)*sin(x)];
%   PointA_med(i+1,:) = [r2_A(i+1)*cos(x), r2_A(i+1)*sin(x)];
%   PointB_int(i+1,:) = [r1_B(i+1)*cos(x), r1_B(i+1)*sin(x)];
%   PointB_med(i+1,:) = [r2_B(i+1)*cos(x), r2_B(i+1)*sin(x)];
% end
% c1=polyfit(r1_A,r1_B,14);
% c2=polyfit(r2_A,r2_B,14);
% 
% Theta = linspace(0,2*pi,251);
% 
% f1 = polyval(c1,Theta);
% f2 = polyval(c2,Theta);
% 
% figure;
% plot(Theta,r1_B)