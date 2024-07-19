clear 
close all

EXP = 4;

%EXP2_AWI.mat  EXP2_HO.mat  EXP4_AWI.mat  EXP4_Dartmouth.mat
nameList = {'AWI (HO)', 'Dart (SSA)', 'Dart (SSA, 2.5km)'};
fileList = {'EXP4_AWI.mat', 'EXP4_Dartmouth.mat', 'EXP4_2_5_Dartmouth'};

Nf = numel(fileList);
% load solutions
for i = 1: Nf
	temp{i} = load(fileList{i});
end

% plot front position
figure('Position', [0, 800, 1000, 400])
for i = 1: Nf
	subplot(1,2,1)
	plot(temp{i}.time, mean(temp{i}.distance(1:4,:)))
	hold on
	subplot(1,2,2)
	plot(temp{i}.time, mean(temp{i}.distance(5:8,:)))
	hold on
end
subplot(1,2,1)
xlim([0,1000])
ylim([200, 500]*1e3)
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Front position (m)')
title('Caprona')
legend(nameList,'location', 'best')

subplot(1,2,2)
xlim([0,1000])
ylim([450, 750]*1e3)
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Front position (m)')
legend(nameList,'location', 'best')
title('Halbrane')

% plot front position
figure('Position', [0, 800, 1000, 400])
for i = 1: Nf
   subplot(1,2,1)
   plot(temp{i}.time, std(temp{i}.distance(1:4,:)))
   hold on
   subplot(1,2,2)
   plot(temp{i}.time, std(temp{i}.distance(5:8,:)))
   hold on
end
subplot(1,2,1)
xlim([0,1000])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Std Front position (m)')
title('Caprona')
legend(nameList,'location', 'best')

subplot(1,2,2)
xlim([0,1000])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Std Front position (m)')
legend(nameList,'location', 'best')
title('Halbrane')
return



%
figure('Position', [0, 800, 1000, 400])
for i = 1: 4
	subplot(1,2,1)
   plot(temp{i}.time, std(temp{i}.distance))
   hold on
end
xlim([0,1000])
xlabel('Time (a)')
set(gcf,'Color','w');
ylabel('Std Front position (m)')
legend(nameList,'location', 'best')

% plot front vel
figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, mean(temp{i}.vel))
   hold on
end
xlim([0,1000])
ylim([0,800])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Frontal Velocity (m/a)')
legend(nameList,'location', 'best')

figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, std(temp{i}.vel))
   hold on
end
xlim([0,1000])
xlabel('Time (a)')
set(gcf,'Color','w');
ylabel('Std Frontal Velocity (m/a)')
legend(nameList,'location', 'best')

% plot front thickness
figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, mean(temp{i}.thickness))
   hold on
end
xlim([0,1000])
ylim([0,350])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Frontal ice thickness (m)')
legend(nameList,'location', 'best')

figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, std(temp{i}.thickness))
   hold on
end
xlim([0,1000])
xlabel('Time (a)')
set(gcf,'Color','w');
ylabel('Std Frontal ice thickness (m)')
legend(nameList,'location', 'best')
