function [a0,a,b] = seg(X,Y,scale_factor)
% X and Y are the x and y data coordinates on a closed boundary
% a0, a and b are the 2N+1 fourier coefficients so that the best-fit
% trig poly has the form
% a0 + a(1)*cos(theta) + a(2)*cos(2*theta) + ...
%          b(1)*sin(theta) + b(2)*sin(2*theta) + ...
%
% The Fourier coefficients are scaled up (or down) by scale_factor.
% Use scale_factor to create larger or smaller versions of the
% same mesh.

N = 8; % 
% plot(X,Y,'r-');
for i=1:length(X)
    r(i) = sqrt(X(i)^2 + Y(i)^2);
    theta(i) = atan2(Y(i),X(i));
end

X0 = [mean(r) ; zeros(2*N,1)];
data = [theta' r'];
options = optimset('MaxFunEvals',300000,'MaxIter',40000,'TolFun',1e-10);
X_s = lsqnonlin(@(X) objective_fun(X,data),X0,[],[],options);


theta_f = linspace(-pi,pi,200);
a0 = X_s(1)*scale_factor;
a = X_s(2:N+1)*scale_factor; 
b = X_s(N+2:2*N+1)*scale_factor;
% 
% hold on
% subplot(1,2,1); plot(theta,r,'r.',theta_f,trig_poly(theta_f,a0,a,b),'k-');
% xlabel('theta'); ylabel('r(theta)');
% XX = trig_poly(theta_f,a0,a,b) .* cos(theta_f);
% YY = trig_poly(theta_f,a0,a,b) .* sin(theta_f);
% subplot(1,2,2); plot(X,Y,'r.',XX,YY,'k-');
% xlabel('x'); ylabel('y');

end

function out = objective_fun(X,data)
% r,theta,a0,a,b
N = (length(X)-1)/2;
a0 = X(1);
a = X(2:N+1);
b = X(N+2:2*N+1);
% a(1) = X(2); a(2) = X(3); a(3) = X(4); a(4) = X(5); a(5) = X(6);
% b(1) = X(7); b(2) = X(8); b(3) = X(9); b(4) = X(10); b(5) = X(11);
r = data(:,2); theta = data(:,1);

out = sum(( trig_poly(theta,a0,a,b) - r ).^2);

end

function out = trig_poly(theta,a0,a,b)

N = length(a);
out = a0;
for i=1:N
    out = out + a(i)*cos(i*theta) + b(i)*sin(i*theta);
end

end

 
