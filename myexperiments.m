clear
close all

today = datestr(date(), 'yyyymmdd');

experiments = [3];
flowmodel = 'MOLHO';

if any(experiments == 1) % exp 1: spin up on a coarse mesh dx=20km {{{
	steps = [1:5];
	resolution = 20e3;
	T = 10000;
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
	relaxT = 2000;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'relaxation time', 2000);
end %}}}
