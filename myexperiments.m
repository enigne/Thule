clear
close all

today = datestr(date(), 'yyyymmdd');

experiments = [10, 12];
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
if any(experiments == 8) % exp 8: relaxation on 2km mesh with maxite=1, discovery {{{
	steps = [8];
	resolution = 2e3;
	relaxT = 10000;
	savePath = [today, '_pseudo_relaxation_2km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 9) % exp 9: project from 2km to 1km mesh and reinitialize {{{
	steps = [1:4, 6];
	coarse_resolution = 2e3;
	resolution = 1e3;
	md = runme('steps', steps, 'resolution', resolution, 'flow model', flowmodel, 'coarse resolution', coarse_resolution);
end %}}}
if any(experiments == 10) % exp 10: relaxation on 1km mesh with maxite=1, discovery {{{
	steps = [8];
	resolution = 1e3;
	relaxT = 2000;
	savePath = [today, '_pseudo_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
if any(experiments == 11) % exp 11: download 1km results from discovery {{{
	steps = [9];
	resolution = 1e3;
	savePath = ['20220803_pseudo_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', 'SSA');
end %}}}
if any(experiments == 12) % exp 12: relaxation on 1km mesh, discovery {{{
	steps = [7];
	resolution = 1e3;
	relaxT = 2000;
	savePath = [today, '_relaxation_1km'];
	md = runme('steps', steps, ...
		'savePath', [savePath],...
		'cluster name', 'discovery',...
		'resolution', resolution, 'flow model', flowmodel, 'relaxation time', relaxT);
end %}}}
