F = load('test_data/Data.mat');

X = F.Data.BrainImaging.data;

O = HeterogeneousID(X,'K',3,'Niter',10000);

plot(O.Z,'o');
ylim([-1 6]);
print('../Results/Z','-dpng');
save('Results.mat','O');
