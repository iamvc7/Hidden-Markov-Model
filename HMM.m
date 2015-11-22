% Vinay 

clear all;
clc;

run data.m;

T = 1000;

[N,N] = size(A);
[N,M] = size(B);

%% Problem 1 - The evaluation problem

alpha = zeros(T,N);
lenobs = length(obs);

for i=1:N
	alpha(1,i) = Pi(i)*B(i,obs(1));
end

for t=1:lenobs-1
	for j=1:N
		sum = 0;
		for i=1:N
			temp = alpha(t,i)*A(i,j);
			sum = sum+temp;
		end
		alpha(t+1,j) = sum * B(j,obs(t+1));
	end
end

lenobs = int16(fix(lenobs));
seq = 0;

for i=1:N
	seq = seq+ alpha(lenobs,i); % seq denotes the probability belongs to the given HMM
end

beta = zeros(lenobs,N);

for i=1:N
	beta(lenobs,i) = 1;
end

for t=lenobs-1:-1:1
	sum = 0;
	for j=1:N
		for i=1:N
			sum=sum+A(j,i)*B(i,obs(t+1))*beta(t+1,i);
		end
		beta(t,j)=sum;   % beta
	end	
end



%% Problem 2 - The uncovering problem

delta = zeros(T,N);  
gamma = zeros(T,N);

for i=1:N
	delta(1,i) = Pi(i) * B(i,obs(1));
	gamma(1,i) = 0;
end

for i=2:lenobs
	for j=1:N
		maxi = -1;
		si = -1;
		for k=1:N
			if(delta(i-1,k)*A(k,j) > maxi)
				maxi = max(maxi,delta(i-1,k)*A(k,j));
				si = k;
			end
		end
		delta(i,j) = maxi*B(j,obs(i));
		gamma(i,j) = si;
	end
end

max_s = -1;
max_val = -1;

for i=1:N
	if(delta(lenobs,i)>max_val)
		max_val = delta(lenobs,i);
		max_s = i;
	end
end

opt_seq = zeros(lenobs,1);
opt_seq(1) = max_s;

for i=lenobs:-1:1
	opt_seq(i)=gamma(i,max_s); %opt_seq represents the correct state sequence for the problem
	max_s=gamma(i,max_s);
end



%% Problem 3 - The Re-estimation problem

temp = zeros(N,N);
gamma = zeros(lenobs,N);
epsilon = zeros(lenobs,N,N);
 
for t=lenobs-1:-1:1
	for i=1:N
 		answer = 0;
 		for j=1:N
 			answer = alpha(t,i)*A(i,j) * B(j,obs(t+1)) * beta(t+1,j);
			answer2 = 0;
			for i_n=1:N
				sum1=0;
				for j_n=1:N
					sum1=sum1+alpha(t,i_n) * A(i_n,j_n) * B(j_n,obs(t+1)) * beta(t+1,j_n);
				end
				answer2 = answer2+sum1;
			end
			answer = answer/answer2;
			epsilon(t,i,j) = answer;
 		end
 	end
end

for t=lenobs-1:-1:1
	for i=1:N
		sum = 0;
		for j=1:N
			sum = sum+ epsilon(t,i,j);
		end
		gamma(t,i) = sum;
	end
end

Pi_bar = zeros(1,N);   % Optimal Pi

for i=1:N
	Pi_bar(i) = gamma(1,i);
end

A_bar = zeros(N,N);		% Optimal A

for i=1:N
	for j=1:N
		num = 0;den=0;
		for t=1:lenobs-1
			num = num + epsilon(t,i,j);
			den = den+gamma(t,i);
		end
		A_bar(i,j) = num/den;
	end
end

B_bar = zeros(N,M);  %Optimal B

for j=1:N
	for k = 1:M
		num = 0; den = 0;
		for t=1:lenobs-1
			if obs(t)==k
				num = num + gamma(t,j);
			end
			den = den+gamma(t,j);
		end
		B_bar(j,k) = num/den;
	end
end

clear answer answer2 den i i_n j j_n k lenobs max_s max_val maxi num si sum sum1 t temp;