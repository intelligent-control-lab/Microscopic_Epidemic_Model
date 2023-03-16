% microscopic epidemic model for Lecture 20 in 16-899 ACRL spring 2020

M = 1000; % Number of agents
K = 100; % Number of time steps
figure(1); hold on;
for CASE = 0:5
switch CASE
    case 0
        u = 1/M;
        u_inf = u;
        delay = 0;
    case 1
        u = 2/M;
        u_inf = u;
        delay = 0;
    case 2
        u = 1/M;
        u_inf = u/10;
        delay = 0;
    case 3
        u = 1/M;
        u_inf = u/10;
        delay = 1;
    case 4
        u = 1/M;
        u_inf = u/10;
        delay = 2;
    case 5
        u = 0.1/M;
        u_inf = u;
        delay = 0;
end

trails = 200;
m = zeros(trails, K+1);

for j = 1:trails
    p = [];
    x = zeros(M, K+1);
    x(1,1) = 1;
    for k = 1:K
        m(j,k) = sum(x(:,k));
        if k <= delay
            p = (1-u)^(m(j,k));
        else
            p = (1-u_inf)^m(j,k-delay)*(1-u)^(m(j,k)-m(j,k-delay));
        end
        for i = 1:M
            if x(i,k) == 1
                x(i, k+1) = 1;
            else
                if rand > p
                    x(i, k+1) = 1;
                end
            end
        end
    end
    m(j,k+1) = sum(x(:,k+1));
end
mean_m = zeros(1, K+1);
max_m = zeros(1, K+1);
min_m = zeros(1, K+1);
for k = 1:K+1
    mean_m(k) = mean(m(:,k));
    max_m(k) = max(m(:,k));
    min_m(k) = min(m(:,k));
end
figure(CASE+10)
image(x.*100)
figure(1); hold on;
times = 0:K;
fill([times';flipud(times')],[max_m'; flipud(min_m')],[CASE/5,0.5,1-CASE/5],'linestyle','none','FaceAlpha', 0.2);
plot(0:K,mean_m,'color',[CASE/5,0.5,1-CASE/5],'LineWidth',2)
end
box on
legend("No Intervention u=0.001", "No Intervention u=0.001",...
    "No Intervention u=0.002", "No Intervention u=0.002",...
    "Immediate Isolation u=0.001, u^*=0.0001", "Immediate Isolation u=0.001, u^*=0.0001",...
    "1-Day Delayed Isolation u=0.001, u^*=0.0001", "1-Day Delayed Isolation u=0.001, u^*=0.0001",...
    "2-Day Delayed Isolation u=0.001, u^*=0.0001", "2-Day Delayed Isolation u=0.001, u^*=0.0001",...
    "Lock-down u^*=0.0001", "Lock-down u^*=0.0001")