% macroscopic SI epidemic model for Lecture 20 in 16-899 ACRL spring 2020

x0 = 1e-3;
N = 100;
x = zeros(1,N+1);

figure(1);clf;hold on;
x(1) = x0;
for beta = 0.2:0.1:1
for i = 1:N
    x(i+1) = x(i) + beta * (1-x(i)) * x(i);
end
plot(0:N, x, 'color', [beta, 0.5, 1-beta])
end
legend("\beta = 0.2", "\beta = 0.3", "\beta = 0.4",...
    "\beta = 0.5", "\beta = 0.6", "\beta = 0.7",...
    "\beta = 0.8", "\beta = 0.9", "\beta = 1.0");
box on;
