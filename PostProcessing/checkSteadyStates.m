clear
close all
addpath('../')
projectsettings;

resolutions = [1e3, 2e3, 5e3, 10e3, 20e3];
figure('Position', [500, 500, 1000, 500])
for i = 1:length(resolutions)
    [time, normDvel, normDH] = convergenceSteadyState(resolutions(i), glacier);
	 subplot(1,2,1)
    semilogy(time(2:end), normDvel,'o')
	 hold on
	 subplot(1,2,2)
    semilogy(time(2:end), normDH,'o')
    hold on
end
subplot(1,2,1)
xlabel('t (years)')
ylabel('max(diff(vel))/diff(time)')
title('velocity')
xlim([0,5000])
subplot(1,2,2)
legend({'1km', '2km', '5km', '10km', '20km'})
xlabel('t (years)')
ylabel('max(diff(H))/diff(time)')
title('ice thickness')
xlim([0,5000])
