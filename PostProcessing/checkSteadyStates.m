clear
close all
addpath('../')
projectsettings;

resolutions = [2e3, 5e3, 10e3, 20e3];
for i = 1:length(resolutions)
    [time, normDvel] = convergenceSteadyState(resolutions(i), glacier);
    semilogy(time(2:end), normDvel,'o')
    hold on
end
legend({'2km', '5km', '10km', '20km'})
xlabel('t (years)')
ylabel('max(diff(vel))/diff(time)')
