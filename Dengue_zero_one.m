load('I.mat');

dat = Idata.I(3000:15:end);

z1test(dat)

% Z1TEST  Gottwald-Melbourne 0-1 test for chaos
% 
% Z1TEST(X) is the result of the 0-1 test applied to the vector X.
% Result is near to 0 for non-chaotic data and near 1 for chaotic data.
% see  http://arxiv.org/pdf/0906.1418v1.  Paul Matthews, July 2009.
function kmedian=z1test(x)
s=size(x); if s(2)==1; x=x'; end;
N=length(x); j=[1:N];
t=[1:round(N/10)];
M=zeros(1,round(N/10));
c=pi/5+rand(1,100)*3*pi/5;  % 100 random c values in [pi/5,4pi/5]
for its=1:100
   p=cumsum(x.*cos(j*c(its)));q=cumsum(x.*sin(j*c(its)));
   for n=1:round(N/10); 
      M(n)=mean( (p(n+1:N)-p(1:N-n)).^2 + (q(n+1:N)-q(1:N-n)).^2 )- ...
           mean(x)^2*(1-cos(n*c(its)))/(1-cos(c(its)));
   end
   kcorr(its)=corr(t',M');
end
%plot(c,kcorr,'*');xlabel('c');ylabel('k'); % useful diagnostic plot

% Two crude attempts to check for oversampling:
if (max(x)-min(x) )/mean(abs(diff(x))) > 10 | ...
       median(kcorr(c<mean(c))) - median(kcorr(c>mean(c))) > 0.5
   disp('Warning: data is probably oversampled.')
   disp('Use coarser sampling or reduce the maximum value of c.')
end
kmedian=median(kcorr);
msg = sprintf('K = %.3f', median(kcorr));

figure;
subplot(1,2,1);
plot(t,M, 'LineWidth',2);xlabel('t');ylabel('M');
annotation('textbox', [0.15 0.7 0.2 0.1], 'String', msg, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white');
subplot(1,2,2);
plot(p,q, 'LineWidth', 2);xlabel('p');ylabel('q');

end
