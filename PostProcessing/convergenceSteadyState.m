function [time, normDvel] = convergenceSteadyState(resolution)
glacier = 'Thule';
projPath = ['/Users/gongcheng/Dartmouth/', glacier];
%% load model {{{
flowmodel = 'SSA';
suffix = ['_', num2str(resolution/1000, '%.0f'), 'km'];
if resolution>=20e3
    stepName = 'Spinup_';
else
    stepName = 'Relaxation_';
end
org=organizer('repository', [projPath, '/Models/'], 'prefix', ['Model_' glacier '_'], 'steps', 0);
md =loadmodel(org, [stepName, flowmodel, suffix]);
disp('Loading model done!'); %}}}
%% compute diff vel {{{
Vel = cell2mat({md.results.TransientSolution(:).Vel});
time = cell2mat({md.results.TransientSolution(:).time});
diffVel = diff(Vel,1,2);
normDvel = sqrt(max(diffVel.^2));
%}}}