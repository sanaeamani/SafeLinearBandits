mu1= 20.8/23;
mu2= 1/23;
mu = [mu1;mu2];
B = [-1,0;0,2];
C = 0.7;

d = 2;
ii = -1:0.0001:1;

xx = [1,1,-1,-1];
yy = [1,-1,-1,1];


L = sqrt(2);
S = 1;
R =1; 
lambda = 1;


J = (mu'*B)';

M = [eye(d);-eye(d)];
b = ones(2*d,1);

xstar = linprog(mu,[mu'*B;M],[C;b])


T = 100000;
lambdaminus = C^2/(d+2)/(S^2)/(norm(B)^2);
delta = 0.01;
bethaa = zeros(1,T);
bound = zeros(1,T);
x = zeros(d,1,T);
rr = zeros(d,1);
A(:,:,1) = lambda*eye(d);
mu_hat = zeros(d,1,T);
SS = polyshape([1000,-1000,-1000,1000],[1000,1000,-1000,-10000]);

%tp = floor((sqrt(d)*2*norm(B)*L*T*(R*sqrt(d*log((1+((T-1)*(L^2)/lambda))/delta))+sqrt(lambda)*S)*sqrt(2)/4/sqrt(lambdaminus)/C)^(2/3))

%tp = floor(tdelta(d,C,B,T,mu,lambdaminus,S,L,R,delta,lambda,xstar));
%tp = floor(T.^(2/3))


tp=0;
Npoints = tp; 
X = randn(Npoints,d); 
X = X./sqrt(sum(X.^2,2)); 
E = nthroot(rand(Npoints,1),d);
X = X.*E;
X = ((C/S)*X)*inv(B');
for (t = 1:tp)
    bethaa(t) = (R*sqrt(d*log((1+((t-1)*(L^2)/lambda))/delta))+sqrt(lambda)*S)^2;
    D = [sqrt(d*bethaa(t))*eye(d);-sqrt(d*bethaa(t))*eye(d)];
    V = (A(:,:,t)^(-0.5))*D'+repmat(mu_hat(:,:,t),1,2*d); %V is n*2n
    poly1 = polyshape(V(1,:),V(2,:));
    SS=intersect(poly1,SS);
    x(:,:,t) = X(t,:)';
    l(t) = mu'*x(:,:,t)+normrnd(0,0.01);
    A(:,:,t+1) = A(:,:,t) + x(:,:,t)*x(:,:,t)';
    rr = rr+l(t)*x(:,:,t);
    mu_hat(:,:,t+1)= inv(A(:,:,t+1))*rr;
    
end

for (t=tp+1:T)
    %bethaa(t) = max (128*n* log(t)*log(((t)^2)/delta) , ((8/3)*log(((t)^2)/delta))^2);
    bethaa(t) = (R*sqrt(d*log((1+((t-1)*(L^2)/lambda))/delta))+sqrt(lambda)*S)^2;
    
    X = [sqrt(d*bethaa(t))*eye(d);-sqrt(d*bethaa(t))*eye(d)];
    V = (A(:,:,t)^(-0.5))*X'+repmat(mu_hat(:,:,t),1,2*d); %V is d*2d
    poly1 = polyshape(V(1,:),V(2,:));
    SS=intersect(poly1,SS);
    V = SS.Vertices;
    [q kk] = size(V);
    h = [b;C*ones(q,1)];
    [U,nr] = con2vert([M;V*B],h);  %U is a p*n each row is a vertex
    
    [m nn] = size(U);
    
    P = zeros(m,2*d);
     for j=1:m
     for i=1:2*d
       P(j,i) = (U(j,:)*V(i,:)');
     end
     end
 oo = find(P==min(min(P)));

 if mod(min(oo),m)==0
     x(:,:,t) = U(m,:)';
 else
 
     x(:,:,t) = U(mod(min(oo),m),:)';

 end
   
 

 l(t) = mu'*x(:,:,t)+normrnd(0,0.01);
 A(:,:,t+1) = A(:,:,t) + x(:,:,t)*x(:,:,t)';
 
     rr = rr+l(t)*x(:,:,t);
 mu_hat(:,:,t+1)= inv(A(:,:,t+1))*rr;
 
 if mod(t,100000)==0 || t==tp+1
 
 figure(8)
 
%subplot(2,1,1) 
%plot(polyshape(V(:,1),V(:,2)))
%title('Confidence Region')
%hold on
%txt = 'true mu';
%plot(mu1,mu2,'r*')
%text(mu1,mu2,txt)

%xlim([-1 1])
%ylim([-1 1])
%grid on
%hold on

%subplot(2,1,2) 
x1 =U(:,1);
y1 = U(:,2);
cx = mean(x1);
cy = mean(y1);
a = atan2(y1 - cy,x1 - cx);
[~, order] = sort(a);
x1 = x1(order);
y1 = y1(order);

plot([x1;x1(1)],[y1;y1(1)],'--')
title('Safe Decision Set, \Delta > 0, t=10000, T^\prime = 10000')
hold on
plot([xx,xx(1)],[yy,yy(1)])


hold on
plot(ii,(-J(1)*ii/J(2)) + (C/J(2)));
hold on
txt = 'x*';


plot(xstar(1),xstar(2),'r*')
text(xstar(1),xstar(2),txt)

xlim([-2.5 2.5])
ylim([-2.5 2.5])
grid on
hold on
pause

end
end

%[r] = mybound(R,delta,lambda,T,L,S,d,tp);

