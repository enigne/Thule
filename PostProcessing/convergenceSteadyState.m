function [time, normDvel, normDH] = convergenceSteadyState(resolution, glacier)
	projPath = ['/totten_1/chenggong/', glacier, '/'];
	%% load model {{{
	flowmodel = 'SSA';
	suffix = ['_', num2str(resolution/1000, '%.0f'), 'km'];
	if resolution>=20e3
		stepName = 'Spinup_';
	elseif resolution<=2e3
		stepName = 'Pseudo_Relaxation_';
	else
		stepName = 'Relaxation_';
	end
	if resolution <=1e3
		org=organizer('repository', [projPath, '/Models/20220803_pseudo_relaxation_1km/'], 'prefix', ['Model_' glacier '_'], 'steps', 0);
	else
		org=organizer('repository', [projPath, '/Models/'], 'prefix', ['Model_' glacier '_'], 'steps', 0);
	end
	md =loadmodel(org, [stepName, flowmodel, suffix]);
	disp('Loading model done!'); %}}}
	%% compute diff vel and H {{{
	Vel = cell2mat({md.results.TransientSolution(:).Vel});
	thickness = cell2mat({md.results.TransientSolution(:).Thickness});
	time = cell2mat({md.results.TransientSolution(:).time});
	diffVel = diff(Vel,1,2)./diff(time);
	diffH = diff(thickness,1,2)./diff(time);
	normDvel = sqrt(max(diffVel.^2));
	normDH = max(abs(diffH));
	%}}}
