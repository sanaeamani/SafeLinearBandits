clc; %close all;
clear all;

ns = [25 50 100];
T = 5000;


max_eig_sum_inv = zeros(length(ns),T);
trace_sum_inv = zeros(length(ns),T);
for n = ns
    
    A = [];
    sum_inv = [];
    
    A(:,:,1) = eye(n);

    % % all reapeating
	x = randn(n,1);
    x = repmat(x,1,T);
    
    % % all reapeating plus small noise
	x = randn(n,1);
    x = repmat(x,1,T);
    x = x + 0.02*randn(n,T);
    
    
    % % random
%     x = rand(n,T);
%     x = -10+20*rand(n,T);

%     % % semirandom 
%     random = 1;
%     x = [];
%     for i=1:ceil(T/n)
%         if i==ceil(T/n), N=T-(ceil(T/n)-1)*n;
%         else N=n;
%         end
%         
%         if random
%             x = [x randn(n,N)];
%             random = 0;
%         else
%             x = [x repmat(randn(n,1),1,N)];
%             random = 1;
%         end
%         
%     end
    
    % % first random then repeating
%     x = randn(n,2*n);
%     x = [x repmat(randn(n,1),1,T-size(x,2))];
%     
    
    size(x)
    sum_inv(:,:,1) = zeros(n,n);
    for i=1:T
        A(:,:,i+1) = A(:,:,i) + x(:,i)*x(:,i)';
        sum_inv(:,:,i+1) = sum_inv(:,:,i) + inv(A(:,:,i));
        max_eig_sum_inv(ns==n,i) = max(eig(sum_inv(:,:,i+1)));
        trace_sum_inv(ns==n,i) = trace(sum_inv(:,:,i+1));
    end

end
    
figure
subplot(1,2,1)
plot(1:T,max_eig_sum_inv)
grid on
xlabel ('T')
ylabel ('eig_1(sum(A_t^{-1}))')
% title(sprintf('n=%d',n))
hold all
str = {};
for n=ns
    str = [str , strcat('n = ' , num2str(n))];
end
legend(str);

subplot(1,2,2)
plot(1:T,trace_sum_inv)
grid on
xlabel ('T')
ylabel ('trace(sum(A_t^{-1}))')
% title(sprintf('n=%d',n))
hold all
str = {};
for n=ns
    str = [str , strcat('n = ' , num2str(n))];
end
legend(str);

 