clear
close all

today = datestr(date(), 'yyyymmdd');

experiments = [6, 7];
flowmodel = 'SSA';

if any(experiments == 1) % exp 1: spin up on a coarse mesh dx=20km {{{
	steps = [1:5];
	resolution = 20e3;
	T = 20000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'final time', T);
end %}}}
if any(experiments == 2) % exp 2: project to 10km mesh and reinitialize {{{
	steps = [1:4, 6];
	resolution = 10e3;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel);
end %}}}
if any(experiments == 3) % exp 3: relaxation on 10km mesh {{{
	steps = [7];
	resolution = 10e3;
	relaxT = 10000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 4) % exp 4: project from 10km to 5km mesh and reinitialize {{{
	steps = [1:4, 6];
	coarse_resolution = 10e3;
	resolution = 5e3;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 5) % exp 5: relaxation on 5km mesh {{{
	steps = [7];
	resolution = 5e3;
	relaxT = 10000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 6) % exp 6: project from 5km to 2km mesh and reinitialize {{{
	steps = [1:4, 6];
	coarse_resolution = 5e3;
	resolution = 2e3;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 7) % exp 7: relaxation on 2km mesh {{{
	steps = [7];
	resolution = 2e3;
	relaxT = 5000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
