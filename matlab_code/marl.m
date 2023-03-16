% MARL in the epidemic model for Lecture 21 in 16-899 ACRL spring 2020

M = 50; % Number of Agents
K = 100; % Time steps
u_list = [0, 1/M, 10/M]; %Low, Medium, High

p = [];
Q = 10 .* ones(2, M+1, length(u_list));
epsilon = 0.1;
alpha = 1;
beta = 1;
delta = 0.5;
figure(2);hold on;
n_ep = 200;
for episode = 1:n_ep
    x = zeros(M, K+1);
    x(1,1) = 1;
    m = zeros(1, K+1);
    m(1) = sum(x(:,1));
    epsilon = 0.5*(n_ep - episode)/n_ep;
for k = 1:K
    u = zeros(1,M);
    % Choosing control
    for i = 1:M
        u(i) = greedy(Q,x(i,k),m(k),epsilon);
    end
    % Dynamics
    for i = 1:M
        if x(i, k) == 1
            x(i, k+1) = 1;
            if u(i)~= 0
                for j = 1:M
                    if x(j, k) == 0
                        p = min(u_list(u(i)), u_list(u(j)));
                        sample = rand;
                        if sample < p
                            x(j, k+1) = 1;
                        end
                    end
                end
            end
        end
    end
    m(k+1) = sum(x(:,k+1));
    % Q-Learning
    for i = 1:M
        %l = alpha * exp(1/(u_list(u(i))-1.1)) + m(k+1)/M;
        %l = x(i,k+1) + alpha * exp(1/(u_list(u(i))-1.1));
        l = x(i,k+1) + alpha * exp(1/(u_list(u(i))-1.1)) + x(i,k)*u(i);
        dQ = (l + delta * min(Q(x(i,k+1)+1, m(k+1)+1, :)) - Q(x(i,k)+1,m(k)+1,u(i)));
        Q(x(i,k)+1,m(k)+1,u(i)) = Q(x(i,k)+1,m(k)+1,u(i)) + beta * dQ;
    end
    % Low-pass filtering
    Qnew = Q(:,:,:);
    Qnew(:,1,:) = 0.95 .* Q(:,1,:) + 0.05 .* Q(:,2,:);
    for i = 2:M-1
        Qnew(:,i,:) = 0.9 .* Q(:,i,:) + 0.05 .* Q(:,i-1,:) + 0.05 .* Q(:,i+1,:);
    end
    Qnew(:,M,:) = 0.95 .* Q(:,M,:) + 0.05 .* Q(:,M-1,:);
    Q = Qnew;
    if m(k+1) == M
        break;
    end
end
plot(0:k, m(1:k+1),'color',[0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep]);
end

box on
xlabel("Time")
ylabel("Infected")
%% Plot
figure(1);clf; hold on;
n = length(u_list)
for i = 1:n
    plot(1:M-1, Q(1,2:M,i), 'color',[i/n,1-i/n,i/n])
    plot(1:M-1, Q(2,2:M,i), '--','color',[i/n,1-i/n,i/n])
end
if n == 2
    legend("healthy, low", "sick, low", "healthy, high", "sick, high")
else
    legend("healthy, low", "sick, low", ...
        "healthy, medium", "sick, medium", ...
        "healthy, high", "sick, high")
end

box on

xlabel("Infected")
ylabel("Q value")
axis([1,49,0,7])
%% Q-Learning
function u = greedy(Q,x,m,epsilon)
sample = rand(1);
nu = size(Q,3);
[~, u] = min(Q(x+1,m+1,:));
if sample < epsilon
    u = ceil(rand(1)*nu);
else
    if length(u) > 1
        sample = rand(1);
        u = u(ceil(sample*length(u)));
    end
end
end