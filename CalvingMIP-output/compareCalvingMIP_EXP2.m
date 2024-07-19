clear 
close all

EXP = 1;

%EXP2_AWI.mat  EXP2_HO.mat  EXP4_AWI.mat  EXP4_Dartmouth.mat
nameList = {'AWI (HO)', 'HO'};
fileList = {'EXP2_AWI.mat', 'EXP2_HO.mat'};

Nf = numel(fileList);
% load solutions
for i = 1: Nf
	temp{i} = load(fileList{i});
end

% plot front position
figure('Position', [0, 800, 500, 400])
for i = 1: Nf
	plot(temp{i}.time, mean(temp{i}.distance))
	hold on
end
xlim([0,1000])
ylim([640, 770]*1e3)
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Front position (m)')
legend(nameList,'location', 'best')

figure('Position', [0, 800, 500, 400])
for i = 1: Nf
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

return

load CalvingMIP-output/cmip-AWI-ex2-G5000-refined-pulse-linear-normal-100.mat
mdAwi = md;
load CalvingMIP-output/CalvingMIP_Exp1_HO_5km.mat
mdDart = md;
mdList = [mdAwi, mdDart];

nameList = {'AWI (HO)', 'HO'};
figure('Position', [0, 800, 500, 400])
for i = 1:numel(nameList)
	md = mdList(i);
	subplot(3,1,1)
	plot([md.results.TransientSolution(:).time], [md.results.TransientSolution(:).GroundedArea])
	hold on
	subplot(3,1,2)
	plot([md.results.TransientSolution(:).time], [md.results.TransientSolution(:).FloatingArea])
	hold on
	subplot(3,1,3)
	plot([md.results.TransientSolution(:).time], [md.results.TransientSolution(:).IceVolume])
	hold on
end
subplot(3,1,1)
xlim([0,1000])
ylabel('Grounded Area')
xlabel('Time (a)')
subplot(3,1,2)
xlim([0,1000])
ylabel('Floating Area')
xlabel('Time (a)')
subplot(3,1,3)
xlim([0,1000])
ylabel('Ice Volume')
xlabel('Time (a)')
legend(nameList,'location', 'best')
set(gcf,'Color','w');


id=8+1;plotmodel(md,'data',md.results.TransientSolution(id).Vel- md.results.TransientSolution(id).CalvingCalvingrate,'layer#all',7,'mask', md.results.TransientSolution(id).MaskIceLevelset<0,'xlim',[-8e5,8e5],'ylim',[-8e5,8e5],'caxis',ca,'title',['AWI(HO), time=',num2str(md.results.TransientSolution(id).time)])
id=800;plotmodel(mdDart,'data',mdDart.results.TransientSolution(id).Vel - mdDart.results.TransientSolution(id).CalvingCalvingrate,'layer#all',10,'mask', mdDart.results.TransientSolution(id).MaskIceLevelset<0,'figure',2,'caxis',ca, 'title',['Dart(HO), time=',num2str(mdDart.results.TransientSolution(id).time)])
