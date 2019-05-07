sn = [1000000,10000000,100000000,1000000000,8000000000];
pn = [1000000,10000000,100000000];
primex = [32, 128, 1024, 5000, 20000];
sharedx = [512, 1024, 4096, 8192, 10000];
threadx = [32, 64, 128, 256, 512, 1024];
chunkx = [1000000, 10000000, 100000000, 1000000000];
streamx = [1, 2, 4, 8, 16];
spn = [1000000, 10000000, 100000000, 1000000000, 10000000000, 100000000000];

% Sequential algorithms comparison, max chunk size
ssoa = [0.794806, 0.814649, 1.218418, 5.965007, 71.066256];
ssoe = [0.802801, 0.820656, 1.625887, 10.696934, 87.443248];
speedups = ssoe./ssoa;

% Parallel algorithms comparison 1 stream, max chunk size
psoe = [0.422912, 3.233792, 30.274240];
psoa = [0.742176, 52.670464, 4496.165039];

% Varying number of primes stored, N = 1,000,000,000
primey = [10.305, 10.337, 10.383, 10.410, 13.755];

% Varying shared memory size, N = 1,000,000,000
sharedy = [1438.07, 819.83, 395.929, 311.554, 319.218];

% Varying number of threads per block, N = 1,000,000,000
thready = [540.27, 385.826, 280.094, 332.483, 526.637, 1019.407];

% Varying chunk size, N = 1,000,000,000
chunky = [489.740, 373.431, 323.717, 280.094];

% Varying num streams, N = 1,000,000,000, chunk = 10,000,000
streamy = [385.01, 307.05, 269.307, 296.370148, 243.232];

% SSOE vs PSOE
sy = [0.008797, 0.072762, 0.489131, 3.987757, 41.216, 439.5988];
py = [0.00003788, 0.00006144, 0.000053248, 0.083682, 1.628, 50.660];
speedupp = sy./py;

psoe1 = [0.03788, 0.06144, 0.053248, 83.682, 1154.987061];
npsoe = [13.612, 147.729, 1632.747, 17357.244141, 75552.976562];
speedup1 = npsoe./psoe1;

figure(1);
hold on;
yyaxis left;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
plot(sn,ssoe,'-bo','LineWidth',1.2);
plot(sn,ssoa,'-ko','LineWidth',1.2);
xlabel('N');
ylabel('Time taken/s');

yyaxis right;
ylabel('Speedup');
plot(sn, speedups, '--rx', 'LineWidth',1.2);
legend({'Sieve of Eratosthenes', 'Sieve of Atkin', 'Speedup'},'Location','northwest');
title('Comparison of sequential algorithms');
hold off;

figure(2);
hold on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
plot(pn,psoe,'-bo','LineWidth',1.2);
plot(pn,psoa,'-ko','LineWidth',1.2);
xlabel('N');
ylabel('Time taken/ms');
legend({'Sieve of Eratosthenes', 'Sieve of Atkin'},'Location','northwest');
title('Comparison of parallel algorithms');
hold off;

figure(3);
plot(primex, primey, '-bo','LineWidth',1.2);
xlabel('M');
ylabel('Time taken/s');
title('Time taken vs M');

figure(4);
set(gca, 'YScale', 'log');
plot(sharedx, sharedy, '-bo','LineWidth',1.2);
xlabel('Size of shared memory');
ylabel('Time taken/ms');
title('Time taken vs size of shared memory');

figure(5);
set(gca, 'YScale', 'log');
plot(threadx, thready, '-bo','LineWidth',1.2);
xlabel('Number of threads in a block');
ylabel('Time taken/ms');
title('Time taken vs number of threads in a block');

figure(6);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(chunkx, chunky, '-bo','LineWidth',1.2);
xlabel('Chunk size');
ylabel('Time taken/ms');
title('Time taken vs chunk size');

figure(7);
set(gca, 'YScale', 'log');
plot(streamx, streamy, '-bo','LineWidth',1.2);
xlabel('Number of streams');
ylabel('Time taken/ms');
title('Time taken vs number of streams');

figure(8);
hold on;
yyaxis left;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
plot(spn,sy,'-bo','LineWidth',1.2);
plot(spn,py,'-ko','LineWidth',1.2);
xlabel('N');
ylabel('Time taken/s');

yyaxis right;
ylabel('Speedup');
set(gca, 'YScale', 'log');
plot(spn, speedupp, '--rx', 'LineWidth',1.2);
legend({'Sequential SOE', 'Parallel SOE', 'Speedup'},'Location','northwest');
title('Parallel vs Sequential SOE');
hold off;

figure(9);
hold on;
yyaxis left;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
plot(sn,psoe1,'-bo','LineWidth',1.2);
plot(sn,npsoe,'-ko','LineWidth',1.2);
xlabel('N');
ylabel('Time taken/ms');

yyaxis right;
ylabel('Speedup');
set(gca, 'YScale', 'log');
plot(sn, speedup1, '--rx', 'LineWidth',1.2);
legend({'Optimized Parallel SOE', 'Naive Parallel SOE', 'Speedup'},'Location','northwest');
title('Naive vs Optimized Parallel SOE');
hold off;
